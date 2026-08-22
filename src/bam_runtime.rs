pub fn bounded_hts_threads(threads: usize) -> usize {
    threads.max(1)
}

#[cfg(test)]
mod tests {
    use super::bounded_hts_threads;

    #[test]
    fn hts_threads_are_never_zero() {
        assert_eq!(bounded_hts_threads(0), 1);
        assert_eq!(bounded_hts_threads(1), 1);
        assert_eq!(bounded_hts_threads(8), 8);
    }
}
