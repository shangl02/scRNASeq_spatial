# How to run jupyterlab
1	Login to HPC server using your Pfizer ID and password through noMachine and open a terminal inside it (It’s recommended to use NoMachine because if you quit NoMachine, the running processes would not be interrupted). <br/>

2	Choose a place to create a folder that will be used to store your container. </br>

3	Download docker container using singularity in HPC.

	singularity pull docker://shl198/sc_ppl:202303. 

(After this command, a singularity image file named sc_ppl_202303.sif is stored in the current directory. This is an container for single cell pipelines, you can use any public docker container that have jupyterlab installed). <br/>

4	Build sandbox using singularity in HPC (this is to make sure you can install new packages in the container in HPC) <br/>

	singularity build –sandbox sc_ppl sc_ppl_2023.sif
(This step will  create a folder named sc_ppl) <br/>

5	In noMachine terminal, login to other compute nodes (you can use lsload in terminal to check what nodes are available in HPC, I list two options here, option one is the standard way to run, option two is more easy but cheat way to do): <br/>
* (Standard way to start an interactive bsub) <br/>
	
        bsub -M 20G -n 1 -q medium -m 'hpccpu200.pfizer.com' -Is bash 

(you can check documentation of bsub to request the number of cores and memory you need, -M sets memory, -n sets number of cores, -q set the queue, -m choose the computing node, this way doesn’t apply to all nodes, you need to select the correct nodes). <br/>
* (you can change to other nodes based on the load)
        
        ssh hpccpu200

6   Start the sandbox <br/>
    
    singularity shell -B /:/media --writable sc_ppl. (This step will access the container)

7   Run jupyterlab with out opening browser <br/>
    
    jupyter-lab --no-browser --port 8080
copy the url shown in the terminal starts with: http://localhost:8080/lab?token= <br/>

8   In your own laptop open a terminal that allows ssh connection (one example is mobaxterm), Run command in terminal <br/>
    
    ssh -L 8080:localhost:8080 lis262@hpccpu101.pfizer.com
9   In the browser of your own laptop, open the url in step 7b.
