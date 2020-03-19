#define NOMINMAX

#include <chrono>
#include <random>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <arrayfire.h>
#include <af/defines.h>
#include <af/seq.h>
#include<cmath>

using namespace std;
using namespace af;

// Skips the Byte Order Mark (BOM) that defines UTF-8 in some text files.
void SkipBOM(std::ifstream& in)
{
    char test[3] = { 0 };
    in.read(test, 3);
    if ((unsigned char)test[0] == 0xEF &&
        (unsigned char)test[1] == 0xBB &&
        (unsigned char)test[2] == 0xBF)
    {
        return;
    }
    in.seekg(0);
}

//Reads in data from csv files into a map to create easy access for later
map<string, int> mapDataFromCSV(string filepath) {

    map< string, int > valMap;
    // Open file:

    ifstream in(filepath);
    assert(in && "file did not open successfully");
    SkipBOM(in);

    // Read the file and load the data:
    string key, record;
    int value;
    string field;
    while (getline(in, record)) {
        istringstream s(record);
        getline(s, key, ',');
        getline(s, field, ',');
        value = stoi(field);

        valMap.insert(pair<string, int>(key, value));
    }

    in.close();

    //// printing map for testing
    //map<string, string>::iterator itr;
    //cout << "\nThe map valMap is : \n";
    //cout << "\tKEY\tELEMENT\n";
    //for (itr = valMap.begin(); itr != valMap.end(); ++itr) {
    //    cout << '\t' << itr->first
    //        << '\t' << itr->second << '\n';
    //}

    return valMap;
}

af::array createEnrolmentMatrix(map<string, int>& studentList, map<string, int>& examList, string filepath) {
    // Open file:
    ifstream in(filepath);
    assert(in && "file did not open successfully");
    SkipBOM(in);

    int row = studentList.size();
    int col = examList.size();

    af::array enrolmentMatrix = af::constant(0, row, col);
    // Read the file and load the data :
    string student, exam, record;
    while (getline(in, record)) {
        istringstream s(record);
        string field;
        getline(s, student, ',');
        getline(s, exam, ',');
        enrolmentMatrix(studentList.at(student), examList.at(exam)) += 1;
    }

    af_sync;
    return enrolmentMatrix;
}

af::array createEnrolmentMatrixForAltData(int numOfExams, int numOfStudents, string filepath) {
    // Open file:
    ifstream in(filepath);
    assert(in && "file did not open successfully");
    SkipBOM(in);

    af::array enrolmentMatrix = af::constant(0, numOfStudents, numOfExams);
    // Read the file and load the data :
    string exam, record;
    int student = 0;
    while (getline(in, record)) {
        istringstream s(record);
        string field;
        while (getline(s, exam, ',')) {
            if (exam != "") {
                enrolmentMatrix(student, stoi(exam) - 1) += 1;
            }
        }
        student++;
    }

    af_sync;
    //af_print(enrolmentMatrix);
    return enrolmentMatrix;
}

void normalizeParticleList(af::array& particleList) {


    af::array negativeValues = particleList < 0;
    af::array reset = af::constant(0, particleList.dims(0), particleList.dims(1));
    replace(particleList, !negativeValues, reset);
    af::array normalizationConstant = sum(particleList, 1);

    replace(particleList, !negativeValues, reset);

    af::array zeroVectorCondition = normalizationConstant == 0;

    af::array resetZeros = af::randu(particleList.dims(0), particleList.dims(1));

    gfor(seq j, particleList.dims(0)) {
        particleList(j, span) = zeroVectorCondition(j) * resetZeros(j, span) + (!zeroVectorCondition)(j) * particleList(j, span);
    }

    normalizationConstant = sum(particleList, 1);
    

    gfor(seq i, particleList.dims(0)) {
        for (int j = 0; j < particleList.dims(1); j++) {
            particleList(i, j) = particleList(i, j) / normalizationConstant(i);
        }
    }
    //int intInitialise = 0;
    //int* max = &intInitialise;

    //unsigned int uInitialise = 1U;
    //unsigned* location = &uInitialise;

    //af::max(max, location, sum(particleList, 1));

    //if (*max > 1) {
    //    cout << "PROBABILITY IS MORE THAN 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    //    af_print(normalizationConstant);
    //    af_print(particleList);
    //}
}

void updateParticles(af::array& particleList, af::array& velocityList) {
    particleList = particleList + velocityList;

    //af::array negativeValues = particleList < 0;
    //af::array reset = af::constant(0, particleList.dims(0), particleList.dims(1));

    //replace(particleList, !negativeValues, reset);

    normalizeParticleList(particleList);
}

void calculateVelocity(af::array& particleList, af::array& particlePbestList, af::array& particleGbest, af::array& velocityList, int numOfParticles, float w, float u1, float u2) {
    random_device rand_dev;
    mt19937 generator(rand_dev());
    uniform_real_distribution<float> rand(0.0, 1);

    af::array tiledParticleGbest = tile(particleGbest, numOfParticles, 1);
    velocityList = velocityList * w + u1 * rand(generator) * (particlePbestList - particleList) + u2 * rand(generator) * (tiledParticleGbest - particleList);
    normalizeParticleList(velocityList); 

}

void validateSchedule(af::array& generatedExamSchedule) {
    af::array validSchedule = sum(generatedExamSchedule, 1);

    af::array validateSchedule = validSchedule != 1;
    validateSchedule = validateSchedule.as(f32);

    if (sum(validateSchedule).scalar<float>() != 0) {
        cout << "INVALID SCHEDULE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TERMINATE NOW";
        af_print(validateSchedule);
    }
}

void sampleParticleDistributions(int numOfParticles, int numOfDimensions, int numOfDistributions, af::array& particleList, af::array& generatedExamSchedule) {
    setSeed(time(NULL));

    int intInitialise = 10;
    int* maxVal = &intInitialise;

    unsigned int uInitialise = 1U;
    unsigned* location = &uInitialise;



    af::max(maxVal, location, sum(particleList, 1));



    generatedExamSchedule = constant(0, numOfParticles * numOfDimensions, numOfDistributions);
    af::array probabilities = af::randu(numOfParticles * numOfDimensions, AF_RANDOM_ENGINE_MERSENNE);

    af::array values = af::constant(0, numOfParticles * numOfDimensions);
    gfor(seq i, numOfParticles * numOfDimensions) {
        for (int j = 0; j < numOfDistributions; j++) {
            probabilities(i) = probabilities(i) - particleList(i, j);
            af::array condition = (probabilities(i) < 0);
            generatedExamSchedule(i, j) = condition.as(f32) * 1;
            probabilities(i) = probabilities(i) + condition.as(f32) * 100;
        }
    }
    generatedExamSchedule = af::moddims(reorder(generatedExamSchedule, 2, 1, 0), numOfDistributions, numOfDimensions, numOfParticles).T();


    af_sync;
}

void calculateFitnessFunction(af::array& generatedExamSchedule, af::array& enrolementMatrix, af::array& fitnessPerParticle, af::array& conflictTimeslots, af::array particleList, int numOfParticles, int numOfExams, int numOfTimeslots) {
    //calculate constraint 1
    af::array averageFitness = af::constant(0, 1, 1, numOfParticles);;
    af::array sampledSchedule = generatedExamSchedule;

    af::array sampledConflictZone;

    for (int i = 0; i < 5; i++) {
        sampleParticleDistributions(numOfParticles, numOfExams, numOfTimeslots, particleList, generatedExamSchedule);
        af::array numberOfClashesPerStudent;
        numberOfClashesPerStudent = matmul(enrolementMatrix, generatedExamSchedule);
        numberOfClashesPerStudent -= 0.5;
        numberOfClashesPerStudent = numberOfClashesPerStudent.as(s16);
        numberOfClashesPerStudent = numberOfClashesPerStudent.as(f32); //minusing 1 across the board and removing negatives

        averageFitness = sum(sum(numberOfClashesPerStudent) * 10);
        af::array condition = averageFitness < fitnessPerParticle;

        replace(fitnessPerParticle, !condition, averageFitness);
        sampledConflictZone = sum(numberOfClashesPerStudent, 0);


        gfor(seq i, numOfParticles) {
            sampledSchedule(span, span, i) = condition(i) * generatedExamSchedule(span, span, i) + (!condition)(i) * sampledSchedule(span, span, i);
            conflictTimeslots(0, span, i) = condition(i) * sampledConflictZone(0, span, i) + (!condition)(i) * conflictTimeslots(0, span, i);
        }

    }
    generatedExamSchedule = sampledSchedule;
    //averageFitness /= 5;
    //fitnessPerParticle = averageFitness;
    af_sync;
}

int updateGbest(af::array& pBestfitnessPerParticle, af::array& particleGbest, af::array generatedExamSchedule, af::array particleList, af::array enrolementMatrix, af::array conflictTimeslots, float* gBestFitnes, int numOfExams, int numOfTimeslots, float learningFactor) {
    int intInitialise = 10;
    int* minValue = &intInitialise;

    unsigned int uInitialise = 1U;
    unsigned* location = &uInitialise;

    random_device rand_dev;
    mt19937 generator(rand_dev());

    af::min(minValue, location, pBestfitnessPerParticle);

    if (*minValue < *gBestFitnes) {
        particleGbest = particleList(seq(*location * numOfExams, (*location + 1) * numOfExams - 1), span);
        //cout << "lower fitness found!!!!!\n";
        *gBestFitnes = *minValue;
        generatedExamSchedule = generatedExamSchedule(span, span, *location);
        af::array dimensionTotal = sum(particleGbest, 1);
        //gfor(seq i, numOfExams) {
        //    for (int j = 0; j < numOfTimeslots; j++) {
        //        af::array condition = (generatedExamSchedule(i, j) == 1);
        //        particleGbest(i, j) = particleGbest(i, j) - (1 - learningFactor) * (!condition) * particleGbest(i, j) + condition * (dimensionTotal(i) - particleGbest(i, j)) * (1 - learningFactor);
        //    }
        //}

        gfor(seq i, numOfExams) {
            for (int j = 0; j < numOfTimeslots; j++) {
                af::array condition = (generatedExamSchedule(i, j) == 1);
                particleGbest(i, j) = (learningFactor) * (!condition) * particleGbest(i, j) + condition * ((dimensionTotal(i) - particleGbest(i, j)) * (1.0 - learningFactor) + particleGbest(i, j));
            }
        }
    }


    //   af::array searchGbest = particleGbest;
    //   float lambda = 1 - pow((*gBestFitnes - 1) / *gBestFitnes, 30000.0/5);
    //   cout << lambda;
    //   af::array localSearchBest = af::constant(0,1);
    //   int localIteration = 0;
    //   localSearchBest(0) = *gBestFitnes;
    //   //local search
    //   if (*gBestFitnes > 15000 ) {
    //       while (localSearchBest(0).scalar<float>() >= *gBestFitnes && localIteration < 10) {
    //           gfor(seq i, numOfExams) {
    //               for (int j = 0; j < numOfTimeslots; j++) {
    //                   searchGbest(i, j) = (1 - lambda) * searchGbest(i, j) + lambda * (4 * searchGbest(i, j) * (1 - searchGbest(i, j)));
    //               }
    //           }
    //           normalizeParticleList(searchGbest);
    //           sampleParticleDistributions(1, numOfExams, numOfTimeslots, searchGbest, generatedExamSchedule);
    //           calculateFitnessFunction(generatedExamSchedule, enrolementMatrix, localSearchBest, conflictTimeslots,particleGbest,1,numOfExams,numOfTimeslots);
    //           localIteration++;
    //           
    //           //cout << "local search best " <<localSearchBest.scalar<float>();
    //           //af_print(searchGbest(1,span));
    //       }
    //   }else {
    //       while (localSearchBest(0).scalar<float>() >= *gBestFitnes && localIteration < 10) {
    //           for (int k = 0; k < numOfExams; k++) {
    //               for (int j = 0; j < numOfTimeslots; j++) {
    //                   normal_distribution<double> distribution(searchGbest(k,j).scalar<float>(), 0.038);
    //                   searchGbest(k, j) = (1 - lambda) * searchGbest(k, j) + lambda * (distribution(generator));
    //               }
    //           }
    //           normalizeParticleList(searchGbest);
    //           sampleParticleDistributions(1, numOfExams, numOfTimeslots, searchGbest, generatedExamSchedule);
    //           calculateFitnessFunction(generatedExamSchedule, enrolementMatrix, localSearchBest, conflictTimeslots,searchGbest,1,numOfExams,numOfTimeslots);
    ///*           cout << "lambda :" << lambda;
    //           af_print(searchGbest);
    //           cout << "iterating\n";*/
    //           localIteration++;
    //       }
    //   }
    //   if (localSearchBest(0).scalar<float>() < *gBestFitnes) {
    //       particleGbest = searchGbest;
    //       *gBestFitnes = localSearchBest(0).scalar<float>();
    //       cout << "found a better value\n";
    //       //af_print(particleGbest);
    //   }
    //   cout << "local iteration :" << localIteration;

    af_sync;
    return *location;
}

void updatePBest(af::array& fitnessPerParticle, af::array& pBestfitnessPerParticle, af::array generatedExamSchedule, af::array& particlePbestList, af::array particleList, int numOfParticles, int numOfExams, int numOfTimeslots, float learningFactor) {

    generatedExamSchedule = moddims(reorder(generatedExamSchedule, 1, 0, 2), numOfTimeslots, numOfParticles * numOfExams).T();

    //af::array dimensionTotal = sum(particleList, 1);
    //af::array condition1 = fitnessPerParticle < pBestfitnessPerParticle;
    //replace(pBestfitnessPerParticle, !condition1, fitnessPerParticle);
    //condition1 = moddims(reorder(tile(condition1, 1, numOfExams), 1, 0, 2), 1, numOfParticles * numOfExams).T();

    //gfor(seq i, numOfParticles * numOfExams) {
    //    for (int j = 0; j < numOfTimeslots; j++) {
    //        af::array condition2 = (generatedExamSchedule(i, j) == 1);
    //        particlePbestList(i, j) = (learningFactor)*condition1(i) * (!condition2) * particleList(i, j)
    //            + condition1(i) * condition2 * ((dimensionTotal(i) - particleList(i, j)) * (1.0 - learningFactor) + particleList(i, j))
    //            + (!condition1)(i) * particlePbestList(i, j);
    //    }
    //}

    af::array dimensionTotal = sum(particlePbestList, 1);
    af::array condition1 = fitnessPerParticle < pBestfitnessPerParticle;
    replace(pBestfitnessPerParticle, !condition1, fitnessPerParticle);
    condition1 = moddims(reorder(tile(condition1, 1, numOfExams), 1, 0, 2), 1, numOfParticles * numOfExams).T();

    gfor(seq i, numOfParticles * numOfExams) {
        for (int j = 0; j < numOfTimeslots; j++) {
            af::array condition2 = (generatedExamSchedule(i, j) == 1);
            particlePbestList(i, j) = particlePbestList(i, j) - (1 - learningFactor) * condition1(i) * (!condition2) * particlePbestList(i, j)
                + condition1(i) * condition2 * (dimensionTotal(i) - particlePbestList(i, j)) * (1 - learningFactor);
        }
    }


    af_sync;
}


void perturbConflictZone(af::array& particleList, af::array& conflictTimeslots, af::array& generatedExamSchedule, af::array particlePbestList, af::array velocityList, float defaultValue) {

    af::array conflictExamsSlots = matmul(conflictTimeslots, generatedExamSchedule.T());

    int intInitialise = 10;
    int* max = &intInitialise;

    unsigned int uInitialise = 1U;
    unsigned* location = &uInitialise;

    //cout << "perturbing!!!!!!!!!!!!!\n\n\n";

    for (int i = 0; i < conflictExamsSlots.dims(2); i++) {
        for (int j = 0; j < 10; j++) {
            af::max(max, location, conflictExamsSlots(span, span, i));
            //af_print(particleList((i * conflictExamsSlots.dims(1)) + *location, span));
            //cout << *location << "\n\n";
            //particleList((i*conflictExamsSlots.dims(1)) + *location , span) = af::constant(defaultValue, 1, conflictTimeslots.dims(1));
            particleList((i * conflictExamsSlots.dims(1)) + *location, span) = af::randu(1, conflictTimeslots.dims(1));
            velocityList((i * conflictExamsSlots.dims(1)) + *location, span) = af::randu(1, conflictTimeslots.dims(1));
            //af_print(particleList((i * conflictExamsSlots.dims(1)) + *location, span));
            //for (int k = 0; k < particleList.dims(1); k++) {
            //    particleList(i * conflictExamsSlots.dims(1) + *location, k) = 4 * particleList(i * conflictExamsSlots.dims(1) + *location, k) * (1 - particleList(i * conflictExamsSlots.dims(1) + *location, k));
            //}

            //af_print(particleList((i * conflictExamsSlots.dims(1)) + *location, span));

            conflictExamsSlots(0, *location, i) = 0;

        }
    }




    //for (int i = 0; i < conflictExamsSlots.dims(2); i++) {
    //    while(*max != 0 ){
    //        af::max(max, location, conflictExamsSlots(span, span, i));
    //        //af_print(particleList((i * conflictExamsSlots.dims(1)) + *location, span));
    //        //cout << *location << "\n\n";
    //        //particleList((i*conflictExamsSlots.dims(1)) + *location , span) = af::constant(defaultValue, 1, conflictTimeslots.dims(1));
    //        particleList((i * conflictExamsSlots.dims(1)) + *location, span) = af::randu(1, conflictTimeslots.dims(1));
    //        //af_print(particleList((i * conflictExamsSlots.dims(1)) + *location, span));
    //        conflictExamsSlots(0, *location, i) = 0;

    //    }
    //    *max = 10;
    //}

    normalizeParticleList(particleList);
}

int main() {
    timer overall = timer::start();
    timer initialization = timer::start();
    af::setBackend(AF_BACKEND_CUDA); //controls whether run on gpu or cpu
    //af::setBackend(AF_BACKEND_CPU);
    info(); //prints deviceinfo

    //collect data
    map<string, int> studentList = mapDataFromCSV("C:\\Users\\lingz\\Documents\\y4sem2\\IE4102\\Cuda Projects\\Cuda Datasets\\Nott\\students.csv");
    map<string, int> examList = mapDataFromCSV("C:\\Users\\lingz\\Documents\\y4sem2\\IE4102\\Cuda Projects\\Cuda Datasets\\Nott\\exams.csv");

    //set stats
    //int numOfStudents = studentList.size(); // nott
    //int numOfExams = examList.size(); // nott
    //int numOfTimeslots = 23; //nott
    //int numOfStudents = 18419; // car 92
    //int numOfExams = 543; //car 92
    //int numOfTimeslots = 32;//car92
    int numOfStudents = 941; // yor 83
    int numOfExams = 181; //yor 83
    int numOfTimeslots = 21; //yor 83
    //int numOfStudents = 30032; // pur 93
    //int numOfExams = 2419; //pur 93
    //int numOfTimeslots = 42; //pur 93
    int numOfParticles = 100;

    float defaultValue = 1.0 / numOfTimeslots;

    double w = 0.5;
    double u1 = 1;
    double u2 = 3;

    int numOfIterations = INT_MAX;
    float learningFactor = 0.5;
    float gBestFitness = INT_MAX;
    int bestParticle;

    //generate matrices
    //af::array enrolementMatrix = createEnrolmentMatrixForAltData(numOfExams, numOfStudents, "C:\\Users\\lingz\\Documents\\y4sem2\\IE4102\\Cuda Projects\\Cuda Datasets\\car\\car-92.csv");
    af::array enrolementMatrix = createEnrolmentMatrixForAltData(numOfExams, numOfStudents, "C:\\Users\\lingz\\Documents\\y4sem2\\IE4102\\Cuda Projects\\Cuda Datasets\\yor\\yor-83.csv");
    //af::array enrolementMatrix = createEnrolmentMatrix(studentList, examList, "C:\\Users\\lingz\\Documents\\y4sem2\\IE4102\\Cuda Projects\\Cuda Datasets\\Nott\\enrolements.csv");
    af::array generatedExamSchedule = af::constant(0, numOfExams, numOfTimeslots, numOfParticles);
    af::array fitnessPerParticle = af::constant(INT_MAX, 1, 1, numOfParticles);
    af::array pBestfitnessPerParticle = af::constant(INT_MAX, 1, 1, numOfParticles);
    af::array conflictTimeslots = af::constant(0, 1, numOfTimeslots, numOfParticles);



    //initialise particles
    af::array particleList = af::randu(numOfParticles * numOfExams, numOfTimeslots);
    //af::array particlePbestList = af::randu( numOfParticles * numOfExams, numOfTimeslots);
    //af::array particleGbest = af::randu( numOfExams, numOfTimeslots);
    af::array velocityList = af::randu( numOfParticles * numOfExams, numOfTimeslots);

    //af::array particleList = af::constant( defaultValue, numOfParticles * numOfExams, numOfTimeslots);
    af::array particlePbestList = af::constant(defaultValue, numOfParticles * numOfExams, numOfTimeslots);
    af::array particleGbest = af::constant(defaultValue, numOfExams, numOfTimeslots);
    //af::array velocityList = af::constant(0, numOfParticles * numOfExams, numOfTimeslots);

    cout << "\n\ninitialization complete, begining algorithm now: ";
    printf("elapsed seconds: %g\n\n", timer::stop(initialization));
    //begin algorithm
    int intInitialise = 0;
    int* maxValue = &intInitialise;

    unsigned int uInitialise = 1U;
    unsigned* location = &uInitialise;

    int stuckCount = 0;
    int lowestgBest = 0;
    int perturbcount = 0;
    timer start1 = timer::start();
    int perturbMultiplier = 0;
    for (int iteration = 0; iteration < numOfIterations; iteration++) {
        
        timer start2 = timer::start();
        updateParticles(particleList, velocityList);

        af::max(maxValue, location, isNaN(particleList));
        if (*maxValue > 0) {
            cout << "nan detected resetting!!!!!! location = " << *location << "\n";
            //af_print(particleList);
            //break;
            
            particleList = af::randu(numOfParticles * numOfExams, numOfTimeslots);
            normalizeParticleList(particleList);
            //particlePbestList = af::constant(defaultValue, numOfParticles * numOfExams, numOfTimeslots);
            //particleGbest = af::constant(defaultValue, numOfExams, numOfTimeslots);
            //velocityList = af::constant(0, numOfParticles * numOfExams, numOfTimeslots);
        }
        calculateFitnessFunction(generatedExamSchedule, enrolementMatrix, fitnessPerParticle, conflictTimeslots,particleList,numOfParticles, numOfExams, numOfTimeslots);
        updatePBest(fitnessPerParticle, pBestfitnessPerParticle, generatedExamSchedule, particlePbestList, particleList, numOfParticles, numOfExams, numOfTimeslots, learningFactor);
        bestParticle = updateGbest(pBestfitnessPerParticle, particleGbest, generatedExamSchedule, particleList, enrolementMatrix, conflictTimeslots, &gBestFitness, numOfExams, numOfTimeslots, learningFactor);


        if (gBestFitness == 0) {
            break;
        }

        if (lowestgBest == gBestFitness) {
            stuckCount++;
            perturbcount++;

            if (stuckCount > 5) {

                //cout << "scattering particles!!!!!!!!!!!\n\n";
 /*               setSeed(time(NULL));
                velocityList = af::randu(numOfParticles * numOfExams, numOfTimeslots);
                normalizeParticleList(velocityList);*/
                //setSeed(time(NULL));
                //particleList = af::randu(numOfParticles * numOfExams, numOfTimeslots);
                //normalizeParticleList(particleList);
            }



        }else {
            stuckCount = 0;
            perturbcount = 0;
            lowestgBest = gBestFitness;
            perturbMultiplier = 0;
        }

        calculateVelocity(particleList, particlePbestList, particleGbest, velocityList, numOfParticles, w, u1, u2);

        if (stuckCount > 400) {
            break;
        }

        

        if (perturbcount > 10 - perturbMultiplier) {
            cout << "perturbing!!!!!!!!!!!!!!!!! \n\n";
            perturbConflictZone(particleList, conflictTimeslots, generatedExamSchedule, particlePbestList,velocityList,defaultValue);
            setSeed(time(NULL));
            velocityList = af::randu(numOfParticles * numOfExams, numOfTimeslots);
            normalizeParticleList(velocityList);
            perturbcount = 0;
            //if (perturbMultiplier < 8) {
            //    perturbMultiplier++;
            //}
        }

        cout << "lowest fitness so far " << gBestFitness << "\n";
        printf("elapsed seconds: %g for round %d\n\n", timer::stop(start2), iteration);
    }

    //af_print(particleGbest);
    //af_print(generatedExamSchedule);

    validateSchedule(generatedExamSchedule);

    cout << "lowest fitness found, " << gBestFitness << " , learning factor,  " << learningFactor << ", w , " << w << ", u1 ," << u1 << ", u2 ," << u2 << ", round, " << "\n";
    printf("elapsed seconds: %g for round \n\n", timer::stop(start1));



    //af_print(generatedExamSchedule(span, span, bestParticle));
    //af_print(particleList(seq(bestParticle * numOfExams, (bestParticle + 1) * numOfExams - 1), span));

    //ofstream outfile;
    //outfile.open("output.txt");
    //int bestfitnesssofar = INT_MAX;
    //for (w = 0.1; w <= 0.9; w += 0.1) {
    //    for (u1 = 1; u1 <= 10; u1 += 2) {
    //        for (u2 = 1; u2 < 10; u2 += 2) {
    //            for (learningFactor = 0.5; learningFactor <= 0.9; learningFactor += 0.1) {

    //                for (int j = 0; j < 3; j++) {

    //                    gBestFitness = INT_MAX;


    //                    af::array fitnessPerParticle = af::constant(INT_MAX, 1, 1, numOfParticles);
    //                    af::array pBestfitnessPerParticle = af::constant(INT_MAX, 1, 1, numOfParticles);
    //                    af::array conflictTimeslots = af::constant(0, 1, numOfTimeslots, numOfParticles);



    //                    //initialise particles
    //                    af::array particleList = af::randu(numOfParticles * numOfExams, numOfTimeslots);
    //                    //af::array particlePbestList = af::randu( numOfParticles * numOfExams, numOfTimeslots);
    //                    //af::array particleGbest = af::randu( numOfExams, numOfTimeslots);
    //                    //af::array velocityList = af::randu( numOfParticles * numOfExams, numOfTimeslots);

    //                    //af::array particleList = af::constant( defaultValue, numOfParticles * numOfExams, numOfTimeslots);
    //                    af::array particlePbestList = af::constant(defaultValue, numOfParticles * numOfExams, numOfTimeslots);
    //                    af::array particleGbest = af::constant(defaultValue, numOfExams, numOfTimeslots);
    //                    af::array velocityList = af::constant(0, numOfParticles * numOfExams, numOfTimeslots);

    //                    cout << "\n\ninitialization complete, begining algorithm now: ";
    //                    printf("elapsed seconds: %g\n\n", timer::stop(initialization));
    //                    //begin algorithm
    //                    int intInitialise = 0;
    //                    int* maxValue = &intInitialise;

    //                    unsigned int uInitialise = 1U;
    //                    unsigned* location = &uInitialise;

    //                    int stuckCount = 0;
    //                    int lowestgBest = 0;
    //                    int perturbcount = 0;
    //                    timer start1 = timer::start();
    //                    for (int iteration = 0; iteration < numOfIterations; iteration++) {

    //                        timer start2 = timer::start();
    //                        updateParticles(particleList, velocityList);

    //                        af::max(maxValue, location, isNaN(particleList));
    //                        if (*maxValue > 0) {
    //                            cout << "nan detected resetting!!!!!! location = " << *location << "\n";
    //                            //af_print(particleList);
    //                            //break;

    //                            particleList = af::randu(numOfParticles * numOfExams, numOfTimeslots);
    //                            normalizeParticleList(particleList);
    //                            //particlePbestList = af::constant(defaultValue, numOfParticles * numOfExams, numOfTimeslots);
    //                            //particleGbest = af::constant(defaultValue, numOfExams, numOfTimeslots);
    //                            //velocityList = af::constant(0, numOfParticles * numOfExams, numOfTimeslots);
    //                        }
    //                        calculateFitnessFunction(generatedExamSchedule, enrolementMatrix, fitnessPerParticle, conflictTimeslots, particleList, numOfParticles, numOfExams, numOfTimeslots);
    //                        updatePBest(fitnessPerParticle, pBestfitnessPerParticle, generatedExamSchedule, particlePbestList, particleList, numOfParticles, numOfExams, numOfTimeslots, learningFactor);
    //                        bestParticle = updateGbest(pBestfitnessPerParticle, particleGbest, generatedExamSchedule, particleList, enrolementMatrix, conflictTimeslots, &gBestFitness, numOfExams, numOfTimeslots, learningFactor);


    //                        if (gBestFitness == 0) {
    //                            break;
    //                        }

    //                        if (lowestgBest == gBestFitness) {
    //                            stuckCount++;
    //                            perturbcount++;

    //                            //if (stuckCount > 5) {
    //                            ////    cout << "scattering particles!!!!!!!!!!!\n\n";
    //                            //    setSeed(time(NULL));
    //                            //    velocityList = af::randu(numOfParticles * numOfExams, numOfTimeslots);
    //                            //    normalizeParticleList(velocityList);
    //                            ////    //setSeed(time(NULL));
    //                            ////    //particleList = af::randu(numOfParticles * numOfExams, numOfTimeslots);
    //                            ////    //normalizeParticleList(particleList);
    //                            //}



    //                        }
    //                        else {
    //                            stuckCount = 0;
    //                            perturbcount = 0;
    //                            lowestgBest = gBestFitness;
    //                        }

    //                        calculateVelocity(particleList, particlePbestList, particleGbest, velocityList, numOfParticles, w, u1, u2);

    //                        if (stuckCount > 50) {
    //                            break;
    //                        }
    //                        if (perturbcount > 10) {

    //                            perturbConflictZone(particleList, conflictTimeslots, generatedExamSchedule, particlePbestList, defaultValue);
    //                            perturbcount = 0;
    //                        }

    //                        // cout << "lowest fitness so far " << gBestFitness << "\n";
    //                        //printf("elapsed seconds: %g for round %d\n\n", timer::stop(start2), iteration);
    //                    }
    //                    if (lowestgBest < bestfitnesssofar) {
    //                        bestfitnesssofar = lowestgBest;
    //                    }
    //                    outfile << "lowest fitness found, " << lowestgBest << " , learning factor,  " << learningFactor << ", w , " << w << ", u1 ," << u1 << ", u2 ," << u2 << ", round, " << "\n";
    //                    printf("elapsed seconds: %g ,lowest ftiness found so far %d\n\n", timer::stop(start1), bestfitnesssofar);
    //                }
    //                cout << "\n\n";
    //            }
    //        }
    //    }
    //}



    //outfile.close();






    printf("overall elapsed seconds: %g\n", timer::stop(overall));
    return 0;
}