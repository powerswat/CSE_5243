How to run
	Prerequisite: Please place all the given feature vector files (TFIDF_fin.mat, 
                    bdyFeatMat.txt, bdyVectLabel.txt, tpcFeatMat.txt, 
                    tpcVectLabel.txt, cntVec.txt, and TFnewIdx.txt) in the 
                    directory where all the .m files are located.
	1. Start Matlab at the directory that has all the .m files in it.
	2. Once you started Matlab, type the command "lab2Driver" at the command prompt 
	   and press return.
		>> lab2Driver
        a. Note: The entire execution may take upto about an hour depdending on 
            the stdlinux server status. The process includes txt file reading, 
            transform data into matrix format, Training and testing on those 
            two models (KNN and NB) with two different feature vector sets 
            (Full and Half) and three different train/test split (70/30, 
            80/20, and 90/10).
    3. If you want to check all the TFIDF information for all the terms in the 
        Reuter dataset as a txt file format, you can run "contextToTxt" code as below:
        >> contextToTxt
        a. It might take around a couple of minutes.
        b. The TFIDF information will be written in TFIDF directory.

How to check the experiment results.
	1. There will be three txt files (accuracy.txt, online_efficiency.txt, and 
        offline_efficiency.txt).
    2. Each row represents the train/test split ratio
        applied to the models. Each column represents the model and the dataset
        used for the experiemnt.
    2. accuracy.txt contains all the accuracy information retrieved from all
        the evaluation settings. 
    3. offline_efficiency.txt contains all the offline efficincy information
        (elasped time for training) retrieved from all the evaluation settings. 
    4. online_efficiency.txt contains all the online efficiency information 
        (elapsed time for testing) retrieved from all the evaluation settings.
    5. After you run convertToTxt.m you will see 19043 txt files each of which
        contains a 5-column table. The number in each file name represents
        the order of Reuter article that has body text.
        The data in each column represents as follows:
        a. The first column represents the features in the particular article.
        b. The second column represents the word frequency of the terms in 
            the article
        c. The third column represents the TF value calculated by the TFIDF 
            formula
        d. The fourth column represents the IDF value calculated by the 
            TFIDF formula
        e. The last column represents the final TFIDF value calculated by 
            the TFIDF formula.
