// ParallelProgrammingLabs.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <mpi.h>
#include "AudioFile.h"
#include <ctime>
#include <stdio.h>
unsigned int start_time; // начальное время
unsigned int end_time; // конечное время
//unsigned int search_time = end_time - start_time; // искомое время


AudioFile<double> inputAudioFile;

using namespace std;

int main(int argc, char **argv)
{
	
	int samplesInOne = 10;
	int procRank, procNum, rc,rootProc=0;
	if (rc = MPI_Init(&argc, &argv))
	{
		cout << "Ошибка запуска, выполнение остановлено " << endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	if (procRank == 0)
	{
		cout << "The number of processes: " << procNum << " my number is " << procRank << endl;

		string audioFileName = "trek";
		inputAudioFile.load(audioFileName + ".wav");
		//for (int samplesInOne1 = 1; samplesInOne1 <= 10; samplesInOne1++)
		//int k = 100;
		//start_time = clock();
		//omp_set_num_threads(k);

		AudioFile<double> audioFile = inputAudioFile;

		int sampleRate = audioFile.getSampleRate();
		int bitDepth = audioFile.getBitDepth();

		int numSamples = audioFile.getNumSamplesPerChannel();
		double lengthInSeconds = audioFile.getLengthInSeconds();

		int numChannels = 1;
		
		AudioFile<double> compressedAudioFile;
		int numSamplesCompressed = numSamples / (samplesInOne);
		compressedAudioFile.setSampleRate(sampleRate / (samplesInOne));
		compressedAudioFile.setBitDepth(bitDepth);
		compressedAudioFile.setNumSamplesPerChannel(numSamplesCompressed);
		compressedAudioFile.setNumChannels(1);
		compressedAudioFile.setAudioBufferSize(numChannels, numSamplesCompressed);

		vector<vector<double>>  newBuffer;
		//double newBuffer1[];
		newBuffer.resize(numChannels);
		for (int i = 0; i < numChannels; i++)
			newBuffer[i].resize(numSamplesCompressed);
		int numSampForProc = numSamples / (procNum - 1) + 1;
		MPI_Bcast(&numSampForProc, 1, MPI_INT, rootProc, MPI_COMM_WORLD);
		MPI_Bcast(newBuffer.data(), newBuffer.size, MPI_DOUBLE, rootProc,MPI_COMM_WORLD);


	}
	else
	{
		for (int i = procRank; i < numSamplesCompressed; i += procNum)
		{


			for (int ii = 0; ii < numChannels; ii++)
			{
				double sum = 0;
				vector<double> compSamples;

				for (int j = 0; j < samplesInOne; j++)
				{
					double tmp = audioFile.samples[ii][i * samplesInOne + j];
					sum += tmp;
					compSamples.push_back(tmp);
				}
				sort(compSamples.begin(), compSamples.end());

				//inputSamplesCount -= intersectionNum;
				sum /= samplesInOne;
				sum += compSamples[compSamples.size() / 2];
				sum /= 2;
				newBuffer[ii][i] = sum;
				//cout << "sum["<< ii <<"][" << i << "]=" << sum << endl;
			}
		}
	}
	

	MPI_Finalize();
	compressedAudioFile.setAudioBuffer(newBuffer);
	compressedAudioFile.save("compressed" + audioFileName + to_string(samplesInOne) + ".wav", AudioFileFormat::Wave);
	//cout<< "runtime for " << nb_threads <<" threads = " << (end_time - start_time) / 1000.0 << endl; // время работы программы  
	cout << "finish" << endl;
	



}
