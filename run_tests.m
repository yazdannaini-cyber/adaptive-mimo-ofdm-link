function run_tests()
%RUN_TESTS Run the complete V15.1 deterministic test suite.
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));
run_all_tests_v15_1();
end
