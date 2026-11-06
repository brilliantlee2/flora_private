use std::collections::VecDeque;

use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use crate::barcode::bounded_levenshtein;

pub fn cluster_directional(
    counts: &HashMap<String, usize>,
    threshold: usize,
) -> HashMap<String, String> {
    if threshold == 0 || counts.len() <= 1 {
        return counts
            .keys()
            .map(|umi| (umi.clone(), umi.clone()))
            .collect();
    }
    if counts.len() <= 64 {
        let graph = adjacency(counts, threshold);
        return create_map_to_correct_umi(connected_components(&graph, counts));
    }

    let umis = counts.keys().collect::<Vec<_>>();
    let abundances = umis.iter().map(|umi| counts[*umi]).collect::<Vec<_>>();
    let mut roots = (0..umis.len()).collect::<Vec<_>>();
    roots.sort_unstable_by(|&a, &b| {
        abundances[b]
            .cmp(&abundances[a])
            .then_with(|| umis[a].cmp(umis[b]))
    });
    let mut trie = ActiveUmiTrie::new(&umis, &abundances);
    let mut seen = vec![false; umis.len()];
    let mut map = HashMap::default();
    let mut queue = VecDeque::new();

    // Traverse the same directed graph implicitly; never materialize its edges.
    for root in roots {
        if seen[root] {
            continue;
        }
        seen[root] = true;
        trie.remove(root);
        queue.push_back(root);
        while let Some(current) = queue.pop_front() {
            map.insert(umis[current].clone(), umis[root].clone());
            // ca >= 2*cb-1 is equivalent to cb <= ceil(ca/2), including cb=0.
            let max_count = abundances[current] / 2 + abundances[current] % 2;
            for next in trie.neighbors(umis[current].as_bytes(), threshold, max_count) {
                seen[next] = true;
                trie.remove(next);
                queue.push_back(next);
            }
        }
    }
    map
}

#[derive(Default)]
struct TrieNode {
    children: Vec<(u8, usize)>,
    parent: usize,
    terminal: Option<usize>,
    active: usize,
    min_count: usize,
}

struct ActiveUmiTrie {
    nodes: Vec<TrieNode>,
    leaves: Vec<usize>,
    abundances: Vec<usize>,
    rows: Vec<Vec<usize>>,
    stack: Vec<(usize, usize)>,
    max_length: usize,
    #[cfg(test)]
    last_query_visits: usize,
}

impl ActiveUmiTrie {
    fn new(umis: &[&String], abundances: &[usize]) -> Self {
        let mut trie = Self {
            nodes: vec![TrieNode {
                min_count: usize::MAX,
                ..TrieNode::default()
            }],
            leaves: Vec::with_capacity(umis.len()),
            abundances: abundances.to_vec(),
            rows: Vec::new(),
            stack: Vec::new(),
            max_length: umis.iter().map(|umi| umi.len()).max().unwrap_or(0),
            #[cfg(test)]
            last_query_visits: 0,
        };
        for (index, umi) in umis.iter().enumerate() {
            let mut node = 0;
            trie.nodes[node].active += 1;
            trie.nodes[node].min_count = trie.nodes[node].min_count.min(abundances[index]);
            for base in umi.bytes() {
                let child = trie.nodes[node]
                    .children
                    .iter()
                    .find_map(|&(label, child)| (label == base).then_some(child));
                node = if let Some(child) = child {
                    child
                } else {
                    let child = trie.nodes.len();
                    trie.nodes.push(TrieNode {
                        parent: node,
                        min_count: usize::MAX,
                        ..TrieNode::default()
                    });
                    trie.nodes[node].children.push((base, child));
                    child
                };
                trie.nodes[node].active += 1;
                trie.nodes[node].min_count = trie.nodes[node].min_count.min(abundances[index]);
            }
            trie.nodes[node].terminal = Some(index);
            trie.leaves.push(node);
        }
        trie
    }

    fn remove(&mut self, index: usize) {
        let mut node = self.leaves[index];
        self.nodes[node].terminal = None;
        loop {
            self.nodes[node].active -= 1;
            if node == 0 {
                break;
            }
            node = self.nodes[node].parent;
        }
    }

    fn neighbors(&mut self, query: &[u8], threshold: usize, max_count: usize) -> Vec<usize> {
        #[cfg(test)]
        {
            self.last_query_visits = 0;
        }
        let mut result = Vec::new();
        if self.nodes[0].active == 0 || self.nodes[0].min_count > max_count {
            return result;
        }
        let threshold = threshold.min(self.max_length.max(query.len()));
        let cap = threshold + 1;
        if self.rows.is_empty() {
            self.rows.push(Vec::new());
        }
        self.rows[0].clear();
        self.rows[0].extend((0..=query.len()).map(|distance| distance.min(cap)));
        if let Some(index) = self.nodes[0].terminal {
            if query.len() <= threshold && self.abundances[index] <= max_count {
                result.push(index);
            }
        }
        self.stack.clear();
        self.stack.push((0, 0));
        while let Some(&(node, child_offset)) = self.stack.last() {
            let Some(&(base, child)) = self.nodes[node].children.get(child_offset) else {
                self.stack.pop();
                continue;
            };
            self.stack.last_mut().unwrap().1 += 1;
            if self.nodes[child].active == 0 || self.nodes[child].min_count > max_count {
                continue;
            }
            #[cfg(test)]
            {
                self.last_query_visits += 1;
            }
            let depth = self.stack.len();
            if self.rows.len() <= depth {
                self.rows.push(Vec::new());
            }
            let (parents, children) = self.rows.split_at_mut(depth);
            let previous = &parents[depth - 1];
            let row = &mut children[0];
            row.resize(query.len() + 1, cap);
            row.fill(cap);
            row[0] = depth.min(cap);
            let mut minimum = row[0];
            let first = depth.saturating_sub(threshold).max(1);
            let last = depth.saturating_add(threshold).min(query.len());
            for column in first..=last {
                row[column] = (row[column - 1] + 1)
                    .min(previous[column] + 1)
                    .min(previous[column - 1] + usize::from(base != query[column - 1]))
                    .min(cap);
                minimum = minimum.min(row[column]);
            }
            if minimum > threshold {
                continue;
            }
            if let Some(index) = self.nodes[child].terminal {
                if row[query.len()] <= threshold && self.abundances[index] <= max_count {
                    result.push(index);
                }
            }
            self.stack.push((child, 0));
        }
        result
    }
}

fn adjacency(counts: &HashMap<String, usize>, threshold: usize) -> HashMap<String, Vec<String>> {
    let mut adj = counts
        .keys()
        .cloned()
        .map(|umi| (umi, Vec::new()))
        .collect::<HashMap<_, _>>();
    let umis = counts.keys().cloned().collect::<Vec<_>>();
    for (i, j) in all_pairs(umis.len()) {
        if i >= j {
            continue;
        }
        let a = &umis[i];
        let b = &umis[j];
        if bounded_levenshtein(a, b, threshold).is_none() {
            continue;
        }
        let ca = counts[a];
        let cb = counts[b];
        if ca >= (cb * 2).saturating_sub(1) {
            if let Some(neighbors) = adj.get_mut(a) {
                neighbors.push(b.clone());
            }
        }
        if cb >= (ca * 2).saturating_sub(1) {
            if let Some(neighbors) = adj.get_mut(b) {
                neighbors.push(a.clone());
            }
        }
    }
    adj
}

#[cfg(test)]
fn candidate_pairs(umis: &[String], threshold: usize) -> Vec<(usize, usize)> {
    if umis.len() < 2 {
        return Vec::new();
    }
    if threshold == 0 || umis.len() <= 64 {
        return all_pairs(umis.len());
    }
    let same_len = umis
        .first()
        .map(|umi| umi.len())
        .map(|len| umis.iter().all(|umi| umi.len() == len))
        .unwrap_or(true);
    if !same_len {
        return all_pairs(umis.len());
    }

    let len = umis[0].len();
    let segments = threshold + 1;
    if len < segments {
        return all_pairs(umis.len());
    }

    let bounds = segment_bounds(len, segments);
    let mut buckets: HashMap<(usize, Vec<u8>), Vec<usize>> = HashMap::default();
    for (umi_idx, umi) in umis.iter().enumerate() {
        let bytes = umi.as_bytes();
        for (segment_idx, (start, end)) in bounds.iter().copied().enumerate() {
            buckets
                .entry((segment_idx, bytes[start..end].to_vec()))
                .or_default()
                .push(umi_idx);
        }
    }

    let mut seen_pairs = HashSet::default();
    let mut pairs = Vec::new();
    // One segment survives threshold edits intact, but indels can shift it.
    for (i, umi) in umis.iter().enumerate() {
        let bytes = umi.as_bytes();
        for (segment_idx, (start, end)) in bounds.iter().copied().enumerate() {
            let segment_len = end - start;
            let first = start.saturating_sub(threshold);
            let last = start.saturating_add(threshold).min(len - segment_len);
            for offset in first..=last {
                let key = (segment_idx, bytes[offset..offset + segment_len].to_vec());
                let Some(indices) = buckets.get(&key) else {
                    continue;
                };
                for &j in indices {
                    if i == j {
                        continue;
                    }
                    let pair = if i < j { (i, j) } else { (j, i) };
                    let packed = ((pair.0 as u64) << 32) | pair.1 as u64;
                    if seen_pairs.insert(packed) {
                        pairs.push(pair);
                    }
                }
            }
        }
    }
    pairs
}

fn all_pairs(n: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::with_capacity(n.saturating_mul(n.saturating_sub(1)) / 2);
    for i in 0..n {
        for j in i + 1..n {
            pairs.push((i, j));
        }
    }
    pairs
}

#[cfg(test)]
fn segment_bounds(len: usize, segments: usize) -> Vec<(usize, usize)> {
    let mut bounds = Vec::with_capacity(segments);
    for idx in 0..segments {
        let start = idx * len / segments;
        let end = (idx + 1) * len / segments;
        bounds.push((start, end));
    }
    bounds
}

fn connected_components(
    graph: &HashMap<String, Vec<String>>,
    counts: &HashMap<String, usize>,
) -> Vec<Vec<String>> {
    let mut nodes = graph.keys().cloned().collect::<Vec<_>>();
    nodes.sort_by(|a, b| counts[b].cmp(&counts[a]).then_with(|| a.cmp(b)));
    let mut seen = HashSet::default();
    let mut components = Vec::new();

    for node in nodes {
        if !seen.insert(node.clone()) {
            continue;
        }
        let mut queue = VecDeque::from([node.clone()]);
        let mut component = vec![node];
        while let Some(curr) = queue.pop_front() {
            if let Some(neighbors) = graph.get(&curr) {
                for next in neighbors {
                    if seen.insert(next.clone()) {
                        queue.push_back(next.clone());
                        component.push(next.clone());
                    }
                }
            }
        }
        component.sort_by(|a, b| counts[b].cmp(&counts[a]).then_with(|| a.cmp(b)));
        components.push(component);
    }
    components
}

fn create_map_to_correct_umi(components: Vec<Vec<String>>) -> HashMap<String, String> {
    let mut map = HashMap::default();
    for component in components {
        if let Some(head) = component.first().cloned() {
            for umi in component {
                map.insert(umi, head.clone());
            }
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exhaustive_map(
        counts: &HashMap<String, usize>,
        threshold: usize,
    ) -> HashMap<String, String> {
        let mut graph: HashMap<String, Vec<String>> = counts
            .keys()
            .cloned()
            .map(|umi| (umi, Vec::new()))
            .collect();
        let umis = counts.keys().collect::<Vec<_>>();
        for (i, j) in all_pairs(umis.len()) {
            let (a, b) = (umis[i], umis[j]);
            if bounded_levenshtein(a, b, threshold).is_some() {
                if counts[a] >= (2 * counts[b]).saturating_sub(1) {
                    graph.get_mut(a).unwrap().push(b.clone());
                }
                if counts[b] >= (2 * counts[a]).saturating_sub(1) {
                    graph.get_mut(b).unwrap().push(a.clone());
                }
            }
        }
        create_map_to_correct_umi(connected_components(&graph, counts))
    }

    fn dna(mut value: usize, length: usize) -> String {
        (0..length)
            .map(|_| {
                let base = b"ACGT"[value % 4] as char;
                value /= 4;
                base
            })
            .collect()
    }

    #[test]
    fn active_trie_prunes_abundance_before_visiting_shared_prefixes() {
        let counts = (0..2048)
            .map(|i| (format!("TGCTG{}", dna(i, 10)), 2))
            .collect::<HashMap<_, _>>();
        let umis = counts.keys().collect::<Vec<_>>();
        let abundances = umis.iter().map(|umi| counts[*umi]).collect::<Vec<_>>();
        let mut trie = ActiveUmiTrie::new(&umis, &abundances);
        assert!(trie.nodes.len() <= 1 + 15 * umis.len());
        assert!(trie.neighbors(umis[0].as_bytes(), 3, 1).is_empty());
        assert_eq!(trie.last_query_visits, 0);
        let map = cluster_directional(&counts, 3);
        assert!(map.iter().all(|(umi, parent)| umi == parent));
    }

    #[test]
    fn directional_maps_match_exhaustive_mixed_lengths_and_zero_counts() {
        let counts = [
            "", "A", "AA", "AC", "ACA", "CACA", "NACA", "N!", "!N", "TGCTG", "tgctg",
        ]
        .into_iter()
        .enumerate()
        .map(|(i, umi)| (umi.to_string(), i % 5))
        .collect();
        for threshold in [0, 1, 2, 3, 8, usize::MAX] {
            assert_eq!(
                cluster_directional(&counts, threshold),
                exhaustive_map(&counts, threshold)
            );
        }
    }

    #[test]
    fn directional_maps_match_exhaustive_randomized_dna() {
        let mut state = 174329u64;
        for iteration in 0..40 {
            let mut counts = HashMap::default();
            for _ in 0..100 {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let length = if iteration % 2 == 0 {
                    10
                } else {
                    3 + (state as usize % 10)
                };
                let prefix = if iteration % 3 == 0 { "TGCTG" } else { "" };
                counts.insert(
                    format!("{prefix}{}", dna((state >> 17) as usize, length)),
                    (state >> 33) as usize % 20,
                );
            }
            for threshold in 0..=3 {
                assert_eq!(
                    cluster_directional(&counts, threshold),
                    exhaustive_map(&counts, threshold),
                    "iteration={iteration}, threshold={threshold}"
                );
            }
        }
    }

    #[test]
    fn active_trie_matches_exhaustive_neighbors_after_removal() {
        let counts = (0..256)
            .map(|i| (dna(i, 4), i % 7))
            .collect::<HashMap<_, _>>();
        let umis = counts.keys().collect::<Vec<_>>();
        let abundances = umis.iter().map(|umi| counts[*umi]).collect::<Vec<_>>();
        let mut trie = ActiveUmiTrie::new(&umis, &abundances);
        for i in (0..umis.len()).step_by(3) {
            trie.remove(i);
        }
        for query in ["", "A", "ACGT", "AACGT", "NACG"] {
            for threshold in 0..=3 {
                let mut actual = trie.neighbors(query.as_bytes(), threshold, 3);
                actual.sort_unstable();
                let expected = umis
                    .iter()
                    .enumerate()
                    .filter_map(|(i, umi)| {
                        (i % 3 != 0
                            && abundances[i] <= 3
                            && bounded_levenshtein(query, umi, threshold).is_some())
                        .then_some(i)
                    })
                    .collect::<Vec<_>>();
                assert_eq!(actual, expected, "query={query}, threshold={threshold}");
            }
        }
    }

    #[test]
    fn directional_maps_match_exhaustive_dense_components() {
        for abundance in [0usize, 1, 2, 7] {
            let mut counts = (0..256)
                .map(|i| (format!("TGCTG{}", dna(i, 6)), abundance))
                .collect::<HashMap<_, _>>();
            counts.insert(String::new(), 0);
            counts.insert("TGCTG".to_string(), 19);
            counts.insert("TGCTGN!".to_string(), 2);
            for threshold in 1..=3 {
                assert_eq!(
                    cluster_directional(&counts, threshold),
                    exhaustive_map(&counts, threshold)
                );
            }
        }
    }

    #[test]
    #[ignore = "synthetic timing; run optimized with --release --ignored --nocapture"]
    fn constant_prefix_benchmark() {
        use std::time::Instant;
        for size in [1000, 2000] {
            for abundance in [1usize, 2] {
                let counts = (0..size)
                    .map(|i| (format!("TGCTG{}", dna(i, 10)), abundance))
                    .collect::<HashMap<_, _>>();
                let umis = counts.keys().cloned().collect::<Vec<_>>();
                let started = Instant::now();
                let mut graph = umis
                    .iter()
                    .cloned()
                    .map(|umi| (umi, Vec::new()))
                    .collect::<HashMap<_, _>>();
                // This is the original indexed implementation, including its pair list.
                let pairs = candidate_pairs(&umis, 3);
                let pair_count = pairs.len();
                for (i, j) in pairs {
                    let (a, b) = (&umis[i], &umis[j]);
                    if bounded_levenshtein(a, b, 3).is_some() {
                        if counts[a] >= (2 * counts[b]).saturating_sub(1) {
                            graph.get_mut(a).unwrap().push(b.clone());
                        }
                        if counts[b] >= (2 * counts[a]).saturating_sub(1) {
                            graph.get_mut(b).unwrap().push(a.clone());
                        }
                    }
                }
                let expected = create_map_to_correct_umi(connected_components(&graph, &counts));
                let baseline_seconds = started.elapsed().as_secs_f64();
                let started = Instant::now();
                let actual = cluster_directional(&counts, 3);
                let new_seconds = started.elapsed().as_secs_f64();
                assert_eq!(actual, expected);
                eprintln!("UMI_BENCH size={size} abundance={abundance} pairs={pair_count} baseline_s={baseline_seconds:.6} new_s={new_seconds:.6} equal=true");
            }
        }
    }

    #[test]
    fn directional_clustering_merges_low_count_neighbors() {
        let counts = [
            ("AAAA".to_string(), 10usize),
            ("AAAT".to_string(), 3usize),
            ("TTTT".to_string(), 2usize),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();
        let map = cluster_directional(&counts, 1);
        assert_eq!(map["AAAA"], "AAAA");
        assert_eq!(map["AAAT"], "AAAA");
        assert_eq!(map["TTTT"], "TTTT");
    }

    #[test]
    fn candidate_pairs_keeps_edit_distance_neighbors_via_segment_buckets() {
        let umis = vec![
            "ACGTACGT".to_string(),
            "ACGTTCGT".to_string(),
            "TTTTTTTT".to_string(),
        ];
        let pairs = candidate_pairs(&umis, 1);
        assert!(pairs.contains(&(0, 1)));
    }

    #[test]
    fn indexed_candidates_keep_shifted_indel_neighbors() {
        let mut umis = vec![
            "ACGTCAGTACGATCGTACGA".to_string(),
            "CGTCAGTACGATCGTACGAA".to_string(),
        ];
        for i in 0..63 {
            let suffix = (0..8)
                .map(|bit| if i & (1 << bit) == 0 { 'A' } else { 'C' })
                .collect::<String>();
            umis.push(format!("TTTTTTTTTTTT{suffix}"));
        }
        assert_eq!(bounded_levenshtein(&umis[0], &umis[1], 2), Some(2));
        for threshold in 2..=3 {
            assert!(candidate_pairs(&umis, threshold).contains(&(0, 1)));
        }
    }

    #[test]
    fn clustering_keeps_indel_correction_across_indexing_boundary() {
        let parent = "ACGTCAGTACGATCGTACGA";
        let child = "CGTCAGTACGATCGTACGAA";
        for distinct_umis in [64, 65] {
            let mut counts = HashMap::default();
            counts.insert(parent.to_string(), 100);
            counts.insert(child.to_string(), 1);
            for i in 0..distinct_umis - 2 {
                let suffix = (0..8)
                    .map(|bit| if i & (1 << bit) == 0 { 'A' } else { 'C' })
                    .collect::<String>();
                counts.insert(format!("TTTTTTTTTTTT{suffix}"), 1);
            }
            assert_eq!(cluster_directional(&counts, 2)[child], parent);
        }
    }

    #[test]
    fn indexed_clustering_rejects_distant_segment_matches() {
        let parent = "AAAAAAAAAAAAAAAAAAAA";
        let distant = "AAAAAAAACCCCCCCCCCCC";
        let mut counts = HashMap::default();
        counts.insert(parent.to_string(), 100);
        counts.insert(distant.to_string(), 1);
        for i in 0..63 {
            let suffix = (0..8)
                .map(|bit| if i & (1 << bit) == 0 { 'A' } else { 'C' })
                .collect::<String>();
            counts.insert(format!("TTTTTTTTTTTT{suffix}"), 1);
        }
        let map = cluster_directional(&counts, 2);
        assert_eq!(map[parent], parent);
        assert_eq!(map[distant], distant);
    }

    #[test]
    fn indexed_clustering_preserves_directional_abundance_filter() {
        let parent = "ACGTCAGTACGATCGTACGA";
        let child = "CGTCAGTACGATCGTACGAA";
        let mut counts = HashMap::default();
        counts.insert(parent.to_string(), 100);
        counts.insert(child.to_string(), 51);
        for i in 0..63 {
            let suffix = (0..8)
                .map(|bit| if i & (1 << bit) == 0 { 'A' } else { 'C' })
                .collect::<String>();
            counts.insert(format!("TTTTTTTTTTTT{suffix}"), 1);
        }
        let map = cluster_directional(&counts, 2);
        assert_eq!(map[parent], parent);
        assert_eq!(map[child], child);
    }

    #[test]
    fn indexed_candidates_match_exhaustive_binary_dna_neighbors() {
        let umis = (0..256)
            .map(|value| {
                (0..8)
                    .map(|bit| if value & (1 << bit) == 0 { 'A' } else { 'C' })
                    .collect::<String>()
            })
            .collect::<Vec<_>>();
        for threshold in 1..=3 {
            let candidates = candidate_pairs(&umis, threshold)
                .into_iter()
                .collect::<HashSet<_>>();
            for (i, j) in all_pairs(umis.len()) {
                if bounded_levenshtein(&umis[i], &umis[j], threshold).is_some() {
                    assert!(
                        candidates.contains(&(i, j)),
                        "missing {i}/{j}, distance <= {threshold}"
                    );
                }
            }
        }
    }
}
