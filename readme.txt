Content aware music score binarization
Reproducible Research
email: telmotpinto@gmail.com; ana.rebelo@gmail.com


The root RR folder has two main folders: code and images.

	Code:
	Here is presented all the code impemented:
	- The RUN.m file is the main script that runs all the methods. Comment any method call you don't want to run;
	- The files in capital letters implement the corresponding methods;
	- The files "run_*.m" files are responsible for calling the methods with the correct parameterization and saving the binarizations to the right folders;
	- The "ME_*.m" files implement Misclassification error (also MOPx and FOPx, for the adaptive methods);

	Images:
	This folder has two subdirectories:
	- grey folder has all the 65 grey scale dataset images;
	- ground_truth has the 10 manual binarizations (ground truths) used to calculate the Misclassification Error;


In addition, two other folders are created:

	images/binarizations/
	- where new binarizations are saved;

	results/
	- where optimum thresholds for the global methods and Misclassification Error for all the methods are saved;