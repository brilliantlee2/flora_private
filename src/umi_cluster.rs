use std::collections::VecDeque;

use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use crate::barcode::bounded_levenshtein;

#[cfg(test)]
thread_local! {
    static PAIR_MATERIALIZATIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

pub fn cluster_directional(
    counts: &HashMap<String, usize>,
    threshold: usize,
) -> HashMap<String, String> {
    // Optimize the measured workload only; long sparse UMIs benefit from the
    // legacy selective index rather than scanning every remaining target.
    if threshold == 3 && counts.keys().all(|umi| umi.len() == 15) {
        return cluster_streamed(counts, threshold);
    }
    let graph = adjacency(counts, threshold);
    create_map_to_correct_umi(connected_components(&graph, counts))
}

fn cluster_streamed(counts: &HashMap<String, usize>, threshold: usize) -> HashMap<String, String> {
    let mut nodes = counts
        .iter()
        .map(|(umi, &count)| (umi, count))
        .collect::<Vec<_>>();
    nodes.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
    let mut remaining = (0..nodes.len()).collect::<Vec<_>>();
    let mut positions = remaining.clone();
    let mut map = HashMap::default();
    let mut queue = VecDeque::new();

    for root in 0..nodes.len() {
        if positions[root] == usize::MAX {
            continue;
        }
        remove_remaining(root, &mut remaining, &mut positions);
        map.insert(nodes[root].0.clone(), nodes[root].0.clone());
        queue.push_back(root);

        // A completed component is forward-closed. A queued node will still be
        // expanded, so edges into any already visited node need not be stored.
        while let Some(current) = queue.pop_front() {
            if remaining.is_empty() {
                return map;
            }
            let distance = PreparedDistance::new(nodes[current].0);
            let mut slot = 0;
            while slot < remaining.len() {
                let target = remaining[slot];
                if nodes[current].1 >= (nodes[target].1 * 2).saturating_sub(1)
                    && distance.within(nodes[target].0, threshold)
                {
                    remove_remaining(target, &mut remaining, &mut positions);
                    map.insert(nodes[target].0.clone(), nodes[root].0.clone());
                    queue.push_back(target);
                } else {
                    slot += 1;
                }
            }
        }
    }
    map
}

fn remove_remaining(node: usize, remaining: &mut Vec<usize>, positions: &mut [usize]) {
    let slot = positions[node];
    remaining.swap_remove(slot);
    if slot < remaining.len() {
        positions[remaining[slot]] = slot;
    }
    positions[node] = usize::MAX;
}

struct PreparedDistance<'a> {
    pattern: &'a str,
    masks: [u64; 256],
}

impl<'a> PreparedDistance<'a> {
    fn new(pattern: &'a str) -> Self {
        let mut masks = [0u64; 256];
        if pattern.len() <= 64 {
            for (position, byte) in pattern.bytes().enumerate() {
                masks[byte as usize] |= 1u64 << position;
            }
        }
        Self { pattern, masks }
    }

    fn within(&self, text: &str, threshold: usize) -> bool {
        if self.pattern.len().abs_diff(text.len()) > threshold {
            return false;
        }
        if self.pattern.is_empty() || self.pattern.len() > 64 {
            return bounded_levenshtein(self.pattern, text, threshold).is_some();
        }
        let mut positive = !0u64;
        let mut negative = 0u64;
        let high = 1u64 << (self.pattern.len() - 1);
        let mut score = self.pattern.len();
        for (position, byte) in text.bytes().enumerate() {
            let equal = self.masks[byte as usize];
            let vertical = equal | negative;
            let horizontal = ((equal & positive).wrapping_add(positive) ^ positive) | equal;
            let mut increase = negative | !(horizontal | positive);
            let mut decrease = positive & horizontal;
            if increase & high != 0 {
                score += 1;
            } else if decrease & high != 0 {
                score -= 1;
            }
            increase = (increase << 1) | 1;
            decrease <<= 1;
            positive = decrease | !(vertical | increase);
            negative = increase & vertical;
            // An intermediate score can fall by one per remaining text byte.
            if score.saturating_sub(text.len() - position - 1) > threshold {
                return false;
            }
        }
        score <= threshold
    }
}

fn adjacency(counts: &HashMap<String, usize>, threshold: usize) -> HashMap<String, Vec<String>> {
    let mut adj = counts
        .keys()
        .cloned()
        .map(|umi| (umi, Vec::new()))
        .collect::<HashMap<_, _>>();
    let umis = counts.keys().cloned().collect::<Vec<_>>();
    for (i, j) in candidate_pairs(&umis, threshold) {
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

fn candidate_pairs(umis: &[String], threshold: usize) -> Vec<(usize, usize)> {
    #[cfg(test)]
    PAIR_MATERIALIZATIONS.with(|calls| calls.set(calls.get() + 1));
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

    fn reference_clustering(
        counts: &HashMap<String, usize>,
        threshold: usize,
    ) -> HashMap<String, String> {
        let graph = adjacency(counts, threshold);
        create_map_to_correct_umi(connected_components(&graph, counts))
    }

    #[test]
    fn prepared_distance_matches_dp_exhaustive_short_binary_sequences() {
        let mut sequences = vec![String::new()];
        for length in 1..=6 {
            for value in 0..1usize << length {
                sequences.push(
                    (0..length)
                        .map(|position| {
                            if value & (1 << position) == 0 {
                                'A'
                            } else {
                                'C'
                            }
                        })
                        .collect(),
                );
            }
        }
        for pattern in &sequences {
            let distance = PreparedDistance::new(pattern);
            for text in &sequences {
                for threshold in 0..=6 {
                    assert_eq!(
                        distance.within(text, threshold),
                        bounded_levenshtein(pattern, text, threshold).is_some(),
                        "{pattern:?}/{text:?} threshold={threshold}"
                    );
                }
            }
        }
    }

    #[test]
    fn prepared_distance_matches_dp_word_boundaries_and_arbitrary_utf8_bytes() {
        let sequences = [
            String::new(),
            "A".into(),
            "NacgtN".into(),
            "éACGT".into(),
            "A".repeat(63),
            "A".repeat(64),
            "A".repeat(65),
            format!("{}C", "A".repeat(63)),
            "ACGT".repeat(16),
            format!("CGT{}A", "ACGT".repeat(15)),
            "T".repeat(80),
        ];
        for pattern in &sequences {
            let distance = PreparedDistance::new(pattern);
            for text in &sequences {
                for threshold in [0, 1, 2, 3, 63, 64, 80] {
                    assert_eq!(
                        distance.within(text, threshold),
                        bounded_levenshtein(pattern, text, threshold).is_some(),
                        "lengths={}/{} threshold={threshold}",
                        pattern.len(),
                        text.len()
                    );
                }
            }
        }
    }

    #[test]
    fn prepared_distance_matches_dp_seeded_random_sequences_and_indels() {
        let mut state = 20261113u64;
        for _ in 0..600 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let length = (state % 65) as usize;
            let pattern = (0..length)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGTNacgt"[((state >> 32) % 9) as usize] as char
                })
                .collect::<String>();
            let mut text = pattern.clone();
            if !text.is_empty() {
                text.remove((state as usize) % text.len());
                let position = ((state >> 16) as usize) % (text.len() + 1);
                text.insert(position, 'T');
            }
            text.push('A');
            for threshold in 0..=4 {
                assert_eq!(
                    PreparedDistance::new(&pattern).within(&text, threshold),
                    bounded_levenshtein(&pattern, &text, threshold).is_some()
                );
                assert_eq!(
                    PreparedDistance::new(&text).within(&pattern, threshold),
                    bounded_levenshtein(&text, &pattern, threshold).is_some()
                );
            }
        }
    }

    #[test]
    fn large_fixed_motif_groups_do_not_materialize_candidate_pairs() {
        let counts = (0..256)
            .map(|value| {
                let random = (0..10)
                    .map(|position| b"ACGT"[(value >> (2 * (position % 4))) & 3] as char)
                    .collect::<String>();
                (format!("TGCTG{random}"), 1)
            })
            .collect::<HashMap<_, _>>();
        let expected = reference_clustering(&counts, 3);
        PAIR_MATERIALIZATIONS.with(|calls| calls.set(0));
        let actual = cluster_directional(&counts, 3);
        assert_eq!(actual, expected);
        PAIR_MATERIALIZATIONS.with(|calls| {
            assert_eq!(
                calls.get(),
                0,
                "normal clustering still materializes the candidate graph"
            )
        });
    }

    #[test]
    fn non_fifteen_base_groups_keep_the_legacy_indexed_path() {
        for (length, threshold) in [
            (14, 3),
            (16, 3),
            (64, 3),
            (65, 3),
            (15, 0),
            (15, 2),
            (15, 4),
        ] {
            let counts = [
                ("A".repeat(length), 10),
                (format!("{}C", "A".repeat(length - 1)), 1),
            ]
            .into_iter()
            .collect::<HashMap<_, _>>();
            let expected = reference_clustering(&counts, threshold);
            PAIR_MATERIALIZATIONS.with(|calls| calls.set(0));
            assert_eq!(cluster_directional(&counts, threshold), expected);
            PAIR_MATERIALIZATIONS.with(|calls| assert_eq!(calls.get(), 1));
        }
    }

    #[test]
    fn mixed_and_empty_umi_groups_select_legacy_fallback_at_public_entry() {
        for other in [String::new(), "A".repeat(14), "A".repeat(16)] {
            let counts = [("A".repeat(15), 10), (other, 1)]
                .into_iter()
                .collect::<HashMap<_, _>>();
            let expected = reference_clustering(&counts, 3);
            PAIR_MATERIALIZATIONS.with(|calls| calls.set(0));
            assert_eq!(cluster_directional(&counts, 3), expected);
            PAIR_MATERIALIZATIONS.with(|calls| assert_eq!(calls.get(), 1));
        }
    }

    #[test]
    fn streaming_matches_graph_for_directional_overlaps_and_lexical_ties() {
        let counts = [
            ("AAAA".to_string(), 10),
            ("AAAT".to_string(), 5),
            ("AATT".to_string(), 3),
            ("TTTT".to_string(), 4),
            ("ATTT".to_string(), 1),
            ("TATT".to_string(), 1),
            ("".to_string(), 0),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();
        for threshold in 0..=4 {
            assert_eq!(
                cluster_streamed(&counts, threshold),
                reference_clustering(&counts, threshold)
            );
        }
    }

    #[test]
    fn streaming_matches_graph_for_mixed_byte_lengths_and_reordered_input() {
        let mut counts = HashMap::default();
        let mut state = 20261113u64;
        for length in [0, 1, 3, 4, 15, 16, 63, 64, 65, 80] {
            for _ in 0..9 {
                let sequence = (0..length)
                    .map(|_| {
                        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                        b"ACGTN"[((state >> 32) % 5) as usize] as char
                    })
                    .collect::<String>();
                counts.insert(sequence, (state % 12) as usize);
            }
        }
        counts.insert("éACGT".to_string(), 4);
        for threshold in 0..=3 {
            let expected = reference_clustering(&counts, threshold);
            assert_eq!(cluster_streamed(&counts, threshold), expected);
            let mut entries = counts.iter().collect::<Vec<_>>();
            entries.reverse();
            let reordered = entries
                .into_iter()
                .map(|(key, &value)| (key.clone(), value))
                .collect::<HashMap<_, _>>();
            assert_eq!(cluster_streamed(&reordered, threshold), expected);
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
