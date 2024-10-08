dir = dict()

# system directories
dir["base"] = os.path.join(workflow.basedir, "..")
dir["env"] = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(workflow.basedir, "..", "scripts")
dir["db"] = os.path.join(workflow.basedir, "databases")

# output directories
try:
    assert(config["args"]["output"]) is not None
    dir["out"] = config["args"]["output"]
except (KeyError, AssertionError):
    dir["out"] = "primer_removal_test.out"

# misc output directories
dir["temp"] = os.path.join(dir["out"], "temp")
dir["results"] = os.path.join(dir["out"], "results")
dir["log"] = os.path.join(dir["out"], "logs")
dir["reports"] = os.path.join(dir["out"], "reports")
dir["bench"] = os.path.join(dir["out"], "benchmarks")
dir["output"] = os.path.join(dir["results"], "output")
dir["stats"] = os.path.join(dir["results"], "stats")
