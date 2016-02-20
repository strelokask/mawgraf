#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cmath>

#include "classifier.h"
#include "EasyBMP.h"
#include "linear.h"
#include "argvparser.h"

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;

using CommandLineProcessing::ArgvParser;

typedef vector<pair<BMP*, int> > TDataSet;
typedef vector<pair<string, int> > TFileList;
typedef vector<pair<vector<float>, int> > TFeatures;
//fills MG and arctan;
void Get_Grad_tan(int N, int M, float **MG, float **arctan, BMP* image){
	# define PI 3.14159265359
    for(int i = 1; i+1 < N; i++){
        arctan[i] = new float[M];
        MG[i] = new float[M];
        for(int j = 1; j+1 < M; j++){
			float Gy = 0, Gx = 0, G;
            float f1 = 0.299*image->GetPixel(i-1, j-1).Red + 0.587*image->GetPixel(i-1, j-1).Green + 0.114*image->GetPixel(i-1, j-1).Blue;
            float f2 = 0.299*image->GetPixel(i-1, j).Red + 0.587*image->GetPixel(i-1, j).Green + 0.114*image->GetPixel(i-1, j).Blue;
            float f3 = 0.299*image->GetPixel(i-1, j+1).Red + 0.587*image->GetPixel(i-1, j+1).Green + 0.114*image->GetPixel(i-1, j+1).Blue;
            float f4 = 0.299*image->GetPixel(i, j-1).Red + 0.587*image->GetPixel(i, j-1).Green + 0.114*image->GetPixel(i, j-1).Blue;
            float f6 = 0.299*image->GetPixel(i, j+1).Red + 0.587*image->GetPixel(i, j+1).Green + 0.114*image->GetPixel(i, j+1).Blue;
            float f7 = 0.299*image->GetPixel(i+1, j-1).Red + 0.587*image->GetPixel(i+1, j-1).Green + 0.114*image->GetPixel(i+1, j-1).Blue;
            float f8 = 0.299*image->GetPixel(i+1, j).Red + 0.587*image->GetPixel(i+1, j).Green + 0.114*image->GetPixel(i+1, j).Blue;
            float f9 = 0.299*image->GetPixel(i+1, j+1).Red + 0.587*image->GetPixel(i+1, j+1).Green + 0.114*image->GetPixel(i+1, j+1).Blue;
            Gx = f9 + f7 + 2*f8 - (f1 + 2*f2 + f3);
            Gy = f9 + 2*f6 + f3 - (f1 + 2*f4 + f7);
            G = sqrt(Gx*Gx + Gy*Gy);//modul of gradient
           	arctan[i][j] = (atan2(Gx , Gy)* 180.0) / PI;
            if(arctan[i][j] < 0) arctan[i][j] += 360;
            MG[i][j] = G;//array if modul of gradients
        }	
    }
}
//div to segments
void Div_cnt(int N,int M, int D, int gradus, vector<float>* first_image_features,float** MG,float** arctan){
	for(int x = 1, cntX = 0; x+1 < N; cntX++){
        int rightX, leftX = x, tempX;
        if((leftX + N/D < N) && (cntX+1 < D)){ rightX = leftX + N/D; tempX = 0;}
        else {rightX = N; tempX = 1;}
        x = rightX;
        for(int y = 1, cntY = 0; y+1 < M; cntY++){
            int rightY, leftY = y, tempY;
            if((leftY + M/D < M) && (cntY+1 < D)){ rightY = leftY + M/D; tempY = 0;}
            else {rightY = M; tempY = 1;}
            y = rightY; 
            vector<float> arr(360/gradus);
            for(int i = leftX;  i+tempX < rightX; i++){
                for(int j = leftY; j+tempY < rightY; j++){
                    arr[trunc(arctan[i][j] / gradus)] += MG[i][j];
                }
            }
            float sum = 0.0;
            for(int i=0; i < (360/gradus); i++){ sum += arr[i]*arr[i]; }
            if (sum < 0.00001) for(int i=0; i < (360/gradus); i++) first_image_features->push_back(arr[i]);
            else for(int i=0; i < (360/gradus); i++) first_image_features->push_back(arr[i] / sqrt(sum));}
    }
}

void Get_Color(int N, int M, int D, BMP* image, vector<float>* first_image_features){
	for(int x = 1, cntX = 0; x+1 < N; cntX++){
        int rightX, leftX = x, tempX;
        if((leftX + N/D < N) && (cntX+1 < D)){ rightX = leftX + N/D; tempX = 0;}
        else {rightX = N; tempX = 1;}
        x = rightX;
        for(int y = 1, cntY = 0; y+1 < M; cntY++){
            int rightY, leftY = y, tempY;
            if((leftY + M/D < M) && (cntY+1 < D)){ rightY = leftY + M/D; tempY = 0;}
            else {rightY = M; tempY = 1;}
            y = rightY;
            float f1=0, f2=0, f3=0;
            int kol_vo = 0;
            for(int i = leftX;  i+tempX < rightX; i++){
                for(int j = leftY; j+tempY < rightY; j++){
                	f1 += image->GetPixel(i, j).Red;
                	f2 += image->GetPixel(i, j).Green;
                	f3 += image->GetPixel(i, j).Blue;
	          		kol_vo++;
	          	}
            }
    		first_image_features->push_back(f1/(255*kol_vo));
    		first_image_features->push_back(f2/(255*kol_vo));
    		first_image_features->push_back(f3/(255*kol_vo));
        }
    }
}
// Load list of files and its labels from 'data_file' and
// stores it in 'file_list'
void LoadFileList(const string& data_file, TFileList* file_list) {
    ifstream stream(data_file.c_str());

    string filename;
    int label;
    
    int char_idx = data_file.size() - 1;
    for (; char_idx >= 0; --char_idx)
        if (data_file[char_idx] == '/' || data_file[char_idx] == '\\')
            break;
    string data_path = data_file.substr(0,char_idx+1);
    
 //   while(!stream.eof() && !stream.fail()) {
  while(stream >> filename >> label){//last read twice
//    	stream >> filename >> label;
		if (filename.size())
            file_list->push_back(make_pair(data_path + filename, label));
    }
    stream.close();
}

// Load images by list of files 'file_list' and store them in 'data_set'
void LoadImages(const TFileList& file_list, TDataSet* data_set) {
    for (size_t img_idx = 0; img_idx < file_list.size(); ++img_idx) {
            // Create image
        BMP* image = new BMP();
            // Read image from file
        image->ReadFromFile(file_list[img_idx].first.c_str());
            // Add image and it's label to dataset
        data_set->push_back(make_pair(image, file_list[img_idx].second));
    }
}

// Save result of prediction to file
void SavePredictions(const TFileList& file_list,
                     const TLabels& labels, 
                     const string& prediction_file) {
        // Check that list of files and list of labels has equal size 
    assert(file_list.size() == labels.size());
        // Open 'prediction_file' for writing
    ofstream stream(prediction_file.c_str());

        // Write file names and labels to stream
    for (size_t image_idx = 0; image_idx < file_list.size(); ++image_idx)
        stream << file_list[image_idx].first << " " << labels[image_idx] << endl;
    stream.close();
}

// Exatract features from dataset.
// You should implement this function by yourself =)
void ExtractFeatures(const TDataSet& data_set, TFeatures* features) {
	for (size_t image_idx = 0; image_idx < data_set.size(); ++image_idx) {
        vector<float> first_image_features;
        BMP* image = data_set[image_idx].first;
        int N = image->TellWidth();
        int M = image->TellHeight();
        float **arctan = new float*[N];
        float **MG= new float*[N];
        Get_Grad_tan(N, M, MG, arctan, image);
//9.3-----------------------------------------------------------------
        Div_cnt(N, M, 8, 45, &first_image_features, MG, arctan);
        Div_cnt(N, M, 16, 45, &first_image_features, MG, arctan);
        Div_cnt(N, M, 32, 45, &first_image_features, MG, arctan);
//9.4-----------------------------------------------------------------
		Get_Color(N, M, 8, image, &first_image_features);        
        features->push_back(make_pair(first_image_features, data_set[image_idx].second));
    	for(int i = 1; i+1 < N; i++){ 
    		delete [] arctan[i];
    		delete [] MG[i];
    	}
    	delete [] MG;
        delete [] arctan;
    }
}

// Clear dataset structure
void ClearDataset(TDataSet* data_set) {
        // Delete all images from dataset
    for (size_t image_idx = 0; image_idx < data_set->size(); ++image_idx)
        delete (*data_set)[image_idx].first;
        // Clear dataset
    data_set->clear();
}

// Train SVM classifier using data from 'data_file' and save trained model
// to 'model_file'
void TrainClassifier(const string& data_file, const string& model_file) {
        // List of image file names and its labels
    TFileList file_list; // vector<pair<string, int> >
        // Structure of images and its labels
    TDataSet data_set; // vector<pair<BMP*, int> >
        // Structure of features of images and its labels
    TFeatures features; // vector<pair<vector<float>, int> >
        // Model which would be trained
    TModel model;
        // Parameters of classifier
    TClassifierParams params;
        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);
        // PLACE YOUR CODE HERE
        // You can change parameters of classifier here
    params.C = 0.01;
    TClassifier classifier(params);
        // Train classifier
    classifier.Train(features, &model);
        // Save model to file
    model.Save(model_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

// Predict data from 'data_file' using model from 'model_file' and
// save predictions to 'prediction_file'
void PredictData(const string& data_file,
                 const string& model_file,
                 const string& prediction_file) {
        // List of image file names and its labels
    TFileList file_list;
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // List of image labels
    TLabels labels;	

        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);

        // Classifier 
    TClassifier classifier = TClassifier(TClassifierParams());
        // Trained model
    TModel model;
        // Load model from file
    model.Load(model_file);
        // Predict images by its features using 'model' and store predictions
        // to 'labels'
    classifier.Predict(features, model, &labels);

        // Save predictions
    SavePredictions(file_list, labels, prediction_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

int main(int argc, char** argv) {
    // Command line options parser
    ArgvParser cmd;
        // Description of program
    cmd.setIntroductoryDescription("Machine graphics course, task 2. CMC MSU, 2014.");
        // Add help option
    cmd.setHelpOption("h", "help", "Print this help message");
        // Add other options
    cmd.defineOption("data_set", "File with dataset",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("model", "Path to file to save or load model",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("predicted_labels", "Path to file to save prediction results",
        ArgvParser::OptionRequiresValue);
    cmd.defineOption("train", "Train classifier");
    cmd.defineOption("predict", "Predict dataset");
        
    // Add options aliases
    cmd.defineOptionAlternative("data_set", "d");
    cmd.defineOptionAlternative("model", "m");
    cmd.defineOptionAlternative("predicted_labels", "l");
    cmd.defineOptionAlternative("train", "t");
    cmd.defineOptionAlternative("predict", "p");

        // Parse options
    int result = cmd.parse(argc, argv);

        // Check for errors or help option
    if (result) {
        cout << cmd.parseErrorDescription(result) << endl;
        return result;
    }

        // Get values 
    string data_file = cmd.optionValue("data_set");
    string model_file = cmd.optionValue("model");
    bool train = cmd.foundOption("train");
    bool predict = cmd.foundOption("predict");

        // If we need to train classifier
    if (train)
        TrainClassifier(data_file, model_file);
        // If we need to predict data
    if (predict) {
            // You must declare file to save images
        if (!cmd.foundOption("predicted_labels")) {
            cerr << "Error! Option --predicted_labels not found!" << endl;
            return 1;
        }
            // File to save predictions
        string prediction_file = cmd.optionValue("predicted_labels");
            // Predict data
        PredictData(data_file, model_file, prediction_file);
    }
}