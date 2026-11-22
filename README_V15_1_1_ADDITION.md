# V15.1.1 Hotfix

This update fixes a multi-transmitter dimension mismatch in the causal
pilot-innovation ICI estimator. The full MIMO LOS prior is sliced per Tx before
entering the single-link temporal tracker. It adds Test 123 and the dedicated
`run_tests_v15_1_1` entry point. See `V15_1_1_HOTFIX.md`.
