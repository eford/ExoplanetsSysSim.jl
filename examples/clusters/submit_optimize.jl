function write_pbs(job_basename, run_number)
    #This function writes a PBS script to be submitted to a computing cluster, which runs 'optimize.jl'

    f = open(job_basename*"_"*string(run_number)*".pbs", "w")
    println(f, "#!/bin/tcsh")
    println(f, "#PBS -A ebf11_a_g_sc_default")
    println(f, "#PBS -l nodes=1:ppn=1")
    println(f, "#PBS -l walltime=48:00:00")
    println(f, "#PBS -l pmem=1gb")
    println(f, "#PBS -j oe")
    println(f, "#PBS -m abe")
    println(f, "#PBS -M myh7@psu.edu")
    println(f, "")
    println(f, "cd \$PBS_O_WORKDIR")
    println(f, "")
    println(f, "/gpfs/group/ebf11/default/julia/bin/julia clusters/optimize.jl "*string(run_number))
    close(f)
end





##### NOTE: must run this script from the same directory in which you want to submit the jobs from!

##### To actually write a number of PBS scripts and submit them:

job_basename = "optimize_job" #base name for the job to be submitted
n_submits = 4 #total number of jobs to submit

for i in 1:n_submits
    write_pbs(job_basename, i)
    job_name = job_basename*"_"*string(i)*".pbs"
    run(`qsub $job_name`) #this line submits the job by running 'qsub' in the command line!
    #println(`qsub $job_name`) #this line is just to see what the command looks like in the command line
    println("Job ", i, " submitted.")
end