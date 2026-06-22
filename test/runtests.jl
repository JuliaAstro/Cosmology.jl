using ParallelTestRunner: runtests, find_tests, parse_args
import Cosmology

const init_code = quote
    # Shared code here
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(Cosmology, args; testsuite, init_code)
