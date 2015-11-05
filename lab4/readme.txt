How to run
	Prerequisite: Please place all the given feature vector files (tpcFeat_fin.mat, 
                    bdyFeat_fin.txt) in the directory where all the .m files are located.

	1. Start Matlab at the directory that has all the .m files in it.

	2. Once you started Matlab, type the command "lab2Driver" at the command prompt 
	   and press return.
		>> lab4Driver

        	a. Note: The program supports iterative run by providing multple elements
			using the parameters as the example below. 
			>> lab4Driver([1 2], [0.01 0.05 0.1])

			The first parameter refers to "the number of mininum neighbor points" 
			and the other parameter refers to "the length of epsilon distance". 
			However, it is suffice to check the correctness of the program without 
			providing any parameters.
			
	3. You can control the number of documents to be sampled by manipulating the code
	   inside the program. If you feel that the program execution taking too long for 
	   any reason, please manipulate the code at "line 7" in "genHalfDataset.m" as below:

	   rand_idx = rand_idx(1:3000); (Default: 3000)
		--> rand_idx = rand_idx(1:[the_number_of_documents_you_want_to_sample]);

	   Although this change will degrade the quality of the analysis results, it will
	   derive the final results much faster.

		a. Note: Each clustering model (ex: DBscan using Euclidean distance) 
			should complete execution in 5 ~ 6 minutes with 3000 sampled documents 
			in normal server condition. Therefore the total program execution 
			time with the default setting can take upto 25 minutes under the 
			normal server condition.

	4. All the analysis results will be saved in "result.txt".
