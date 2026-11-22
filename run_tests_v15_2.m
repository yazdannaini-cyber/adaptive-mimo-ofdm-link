function run_tests_v15_2()
%RUN_TESTS_V15_2 Entry point for V15.2 statistical and packet-coding tests.
projectRoot=fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));
run_all_tests_v15_2();
end
