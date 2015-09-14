How to run
	Prerequisite: Please put all the reuter's sgm files in "reuters" directory and
			place the directory where all the .m files are located.
	1. Start Matlab at the directory that has all the .m files in it.
	2. Once you started Matlab, type the command "textPreProc" at the command prompt 
	   and press return.
		>> textPreProc
	3. Once you finished execution you can check the result feature vector by
	   download "FVect.mat" file to your PC or load this using the matlab program by
	   typing the command as below:
		>> load('FVect_fin.mat');
		>> fVec([start_row_num]:[end_row_num], [start_col_num]:[end_col_num])
	   You can display the whole matrix by typing "FVect" at the command prompt but
	   it is not recommenderable since the matrix is too big to be completely 
	   presented in one screen.
	4. You can also check the stemming dictionary generation process by removing 
	   "stemDictionary_fin.mat" file. However, it is also not recommendable since it
	   takes too much time to make due to exploiting online dictionary website.
	5. You can check the TF-IDF caculation process by removing "TFIDF_fin.mat" file.
	6. Stemmed result is saved in bodyTxt_fin.mat. So you can remove this to check
	   how the words are stemmed.
	7. You can check the feature matrix generation process by removing "bdyFeat_fin.mat"
	   "tpcFeat_fin.mat" and "plcFeat_fin.mat"

How to check the feature result.
	1. There are three types of feature vectors generated after the execution.
	   Each row represence the document ID and each column represents the selected
	   features from each category (<topic>, <body>, <places>)
	2. When load "bdyFeat_fin.mat", "bdyVectLabel" represents all the selected features
	   from body texts and "bdyFeatMat" represents all the frequency results of the selected 
	   body feature words from each article.
	3. When load "tpcFeat_fin.mat", "tpcVectLabel" represents all the selected features
	   from topic words and "tpcFeatMat" represents all the frequency results of the selected 
	   topic feature words from each article.
	4. When load "plcFeat_fin.mat", "plcVectLabel" represents all the selected features
	   from place words.
