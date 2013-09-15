//To Compile:
//gcc -o SimpleDSP SimpleDSP.c -lportaudio -lpthread -lm -lrt -Wall -Wextra

/*
*====================Documentation=======================
*					  Future Work
*-Robust localization technique workable around environmental noise
*-Time series filtering (comb filter) [Information Encoding via Frequency]
*-Windowing option (Hann, Hamming etc)
*-Implement additional FFT option (with variable FFT size)
*-Sound Signature Identification (Hybrid Technology Involve Implementation Similar to Speech Recognition)
*-Include calibration function to accurately estimate the frequency range of tone emission (record reference ambient, reference tone, divides and search for narrowest range(while increasing FFT bin size))
*-
*-Time stamp data file
*-Plot 3D graph using gnuplot (RT plotting) [Adjust time axis by v=(max - current), f.ex if v < 3 seconds then advance axis by 3 ]
*-Employ data buffering with dynamic sizing managed by child, 
*-Use shared variable to dynamically control plotting process (1 for entry number, 1 for print (busy looping))
*-Print Procedure (Parent -> Send Data via pipe -> Child recieves and creates data file(stamped by time and date) -> Child thread creates/updates graph)
*-
*-
*=======================Note=============================
*-"Packet" in variable refers to the individual data structure that contains DSP computation results for each individual frequency
*-
*-
*-
*-
*==================End Documentation=====================
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <sys/mman.h>
#include "portaudio.h"
#include "kbhit.h"

#define NANO_SECOND 1000000000
#define LINE_SIZE 128
#define MAX_LOAD 64 //Defines Maximum Number Of Frequencies To Be Plotted In RT

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if 1
#define PA_SAMPLE_TYPE  paFloat32
#define PRINTF_S_TYPE "Float32"
typedef float SAMPLE;
#define SAMPLE_SILENCE  (0.0f)
#define PRINTF_S_FORMAT "%.8f"
#define SAMPLE_TYPE_INDEX 0
#elif 1
#define PA_SAMPLE_TYPE  paInt16
#define PRINTF_S_TYPE "Int16"
typedef short SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#define SAMPLE_TYPE_INDEX 1
#elif 0
#define PA_SAMPLE_TYPE  paInt8
#define PRINTF_S_TYPE "Int8"
typedef char SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#define SAMPLE_TYPE_INDEX 2
#else
#endif

//function protoypes
int compare(const void*, const void*);
int streamCallback(const void*, void*,unsigned long, const PaStreamCallbackTimeInfo*, PaStreamCallbackFlags, void*);
void* DFTFloat32SingleFrequency(void*);
int queueInputDevice();
void buildFrequency();
void destroyFrequency();
int configureInputParameters();
int selectInput();
int parseFrequency();
void DSPdataOutput();
void setupPlot();

typedef struct{
	float* inputVolt;
	int* sampleSize;
	double* currentTime;
} DFTFloat32;

typedef struct{
	float* singlePhase;
	int* bufferSize;
	double* currentTime;
} activeFloat32Mono;

typedef struct{
	float* leftPhase;
	float* rightPhase;
	int* bufferSize;
	double* currentTime;
} activeFloat32Stereo;

typedef struct{
	float* inputSignal;
	double* realF;
	double* imagF;
	double* phaseAngle;
	double* magnitude;
	double* decibel;
	int zeroIndex;
	int frequency;
	int length;
	int samplingRate;
} Float32DFTSingle;

typedef struct {
	int frequency;
	double real;
	double imag;
	double angle;
	double magnitude;
	double decibel;
} DSPPacket;

typedef struct {
	pthread_cond_t cond_IPC;
	pthread_condattr_t condattr_IPC;
	pthread_mutex_t mutex_IPC;
	pthread_mutexattr_t mutexattr_IPC;
    int state;
    int setup;
    double packetTime;
    int packetCount;
    DSPPacket DSPBundle[MAX_LOAD]; //IPC Variable Declaration, <- Non-flexible Implementation...Flexible Implementation -> MAX_LOAD = MAX(SAMPLE_SIZE/2,(*Total Number of Frequencies to be Analayzed*))
    int interval;
    int cond_signal;
    int bufferSize;
} IPCSyncObj;

//default parameters
int SAMPLE_RATE = 44100;
int FRAMES_PER_BUFFER = 512;
int NUM_CHANNELS = 2;

//declare variables
int* deviceIndex;
int* freqList; //stores frequencies to be analyzed
int* freqListSize; //stores the size of the list of frequencies to be analyzed
int* freqIndex; //stores the index of first free element in freqList
int numberOfInput, formatType;
Float32DFTSingle** fq; //double pointer to store a list of Float32DFTSingle* pointers which will be used in the thread workers as argument
PaError paErr;
PaStream* audioStream;
struct timespec* startTime; //define stream start time struct pointer
struct timespec* streamTime; //define active stream time struct pointer
PaStreamParameters inputParameters, outputParameters;

activeFloat32Mono* streamFloat32Mono;
activeFloat32Stereo* streamFloat32Stereo;
DFTFloat32* audioFloat32;

//declare IPC related variables
IPCSyncObj* gnuPlotIPC;

//declare pipe file descriptors
int pipeFD[2]; //file descriptors for pipe
int stdFD[2]; //file descriptors for stdin & stdout

//data output format
int outputFormat = 1; //atomically operated variable to manage data output format

//declare parent thread related variables
pthread_t* t_id;
pthread_attr_t t_attr;
pthread_cond_t cond_buffer;
pthread_condattr_t condattr_buffer;
pthread_mutex_t mutex_buffer;
pthread_mutexattr_t mutexattr_buffer;
//declare child thread related variables
pthread_t child_t_id;
pthread_attr_t child_t_attr;
pthread_cond_t child_cond_buffer;
pthread_condattr_t child_condattr_buffer;
pthread_mutex_t child_mutex_buffer;
pthread_mutexattr_t child_mutexattr_buffer;

//file variables
FILE* fp;

int compare(const void * a, const void * b){
   return (*(int*)a-*(int*)b);
}

int streamCallback(const void *inputBuffer, void *outputBuffer,unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void *userData){
	long unsigned int i;
	switch (formatType){
		case 0:{
			//case paFloat32 sample type
			float* inputStream = (float*)inputBuffer;
			if (NUM_CHANNELS == 1){	//case mono
				activeFloat32Mono* dataBuffer = (activeFloat32Mono*) userData;
				float* phasePtr = (float*)dataBuffer->singlePhase;
				//begin critical section
				pthread_mutex_lock(&mutex_buffer); //acquire lock
				for (i=0; i< framesPerBuffer; i++){
					//printf("singleStream: "PRINTF_S_FORMAT"\n", *phasePtr);
					*phasePtr++ = *inputStream++;
				}
				//timestamp stream audio data
				clock_gettime(CLOCK_MONOTONIC_RAW, streamTime);
				*(dataBuffer->currentTime) = (streamTime->tv_sec - startTime->tv_sec) + ((double)(streamTime->tv_nsec - startTime->tv_nsec)/NANO_SECOND);
				pthread_cond_signal(&cond_buffer);
				//end critical section
				pthread_mutex_unlock(&mutex_buffer); //release lock
			} else if (NUM_CHANNELS == 2){	//case stereo
				activeFloat32Stereo* dataBuffer = (activeFloat32Stereo*) userData;
				float* leftPhasePtr = (float*)dataBuffer->leftPhase;
				float* rightPhasePtr = (float*)dataBuffer->rightPhase;
				//begin critical section
				pthread_mutex_lock(&mutex_buffer); //acquire lock
				for (i=0; i< framesPerBuffer; i++){
					//printf("leftStream: "PRINTF_S_FORMAT" rightStream: "PRINTF_S_FORMAT"\n", *leftPhasePtr, *rightPhasePtr);
					*leftPhasePtr++ = *inputStream++;
					*rightPhasePtr++ = *inputStream++;
				}
				//timestamp stream audio data
				clock_gettime(CLOCK_MONOTONIC_RAW, streamTime);
				*(dataBuffer->currentTime) = (streamTime->tv_sec - startTime->tv_sec) + ((double)(streamTime->tv_nsec - startTime->tv_nsec)/NANO_SECOND);
				pthread_cond_signal(&cond_buffer);
				//end critical section
				pthread_mutex_unlock(&mutex_buffer); //release lock
			}
			break;
		}
		default:
			break;
	}
	return paContinue;
}

//Queue for Audio Input Device and Re-Index for Simplicity
int queueInputDevice(){
	int i;
	int j = 0;
	int inputDeviceCount = 0;
	PaDeviceIndex deviceCount;
	PaDeviceInfo* deviceInfo;
	deviceCount = Pa_GetDeviceCount();
	if (deviceCount < 0) {
		printf("PortAudio Error: %s\n", Pa_GetErrorText((PaError) deviceCount));
		return -1;
	} else if (deviceCount == 0) {
		printf("No Input Devices Available\n");
		return -1;
	} else {
		for (i = 0; i < deviceCount; i++){
			deviceInfo = (PaDeviceInfo*)Pa_GetDeviceInfo(i);
			if (deviceInfo->maxInputChannels > 0) {
				inputDeviceCount++;
			}
		}
		if (inputDeviceCount > 0) {
			numberOfInput = inputDeviceCount;
			deviceIndex = (int*)malloc((numberOfInput) * sizeof(int));
			for (i = 0; i < deviceCount; i++){
				deviceInfo = (PaDeviceInfo*)Pa_GetDeviceInfo(i);
				if (deviceInfo->maxInputChannels > 0) {
					*(deviceIndex + j)= i;
					j++;
				}
			}
		} else {
			printf("No Input Devices Available\n");
			return -1;
		}
	}
	printf("Device Index\tDevice Name\n");
	for (i = 0; i< numberOfInput; i++){
		deviceInfo = (PaDeviceInfo*)Pa_GetDeviceInfo(*(deviceIndex + i));
		printf("%d\t\t%s\n",i, deviceInfo->name);
	}
	return 0;
}

//Configure input parameters
int configureInputParameters(){
	printf("Number of Channels: ");
	scanf("%d", &NUM_CHANNELS);
	printf("Sampling Rate: ");
	scanf("%d", &SAMPLE_RATE);
	printf("Buffer Size: ");
	scanf("%d", &FRAMES_PER_BUFFER);
	return 0;
}

//Configure and test compatibility by using default parameters
int selectInput(){
	PaDeviceIndex index;
	printf("Choose Input Device: ");
	scanf("%d", &index);
	if (index >= numberOfInput) {
		printf("Invalid Index\n");
		return -1;
	}
	printf("Input Device Selected: %s\n", Pa_GetDeviceInfo(*(deviceIndex + index))->name);
	configureInputParameters();
	inputParameters.device = *(deviceIndex + index);
	inputParameters.channelCount = NUM_CHANNELS;
	inputParameters.sampleFormat = PA_SAMPLE_TYPE;
	inputParameters.suggestedLatency = Pa_GetDeviceInfo(inputParameters.device)->defaultLowInputLatency ;
	inputParameters.hostApiSpecificStreamInfo = NULL;
	paErr = Pa_IsFormatSupported(&inputParameters, NULL, SAMPLE_RATE);
	printf("CHANNELS: %d SAMPLE RATE: %d BUFFER SIZE: %d SAMPLE FORMAT: %s\n", NUM_CHANNELS, SAMPLE_RATE, FRAMES_PER_BUFFER, PRINTF_S_TYPE);
	if (paErr != paFormatIsSupported){
		printf("FORMAT REJECTED\n");
		printf("PortAudio Error: %s\n", Pa_GetErrorText(paErr));
		return -1;
	}
	formatType = SAMPLE_TYPE_INDEX;
	printf("FORMAT ACCEPTED\n");
	return 0;
}

//parses the user input and store into a list of frequencies for DFT
int parseFrequency(){
	int i=0, k= 0;
	int c;
	int j, u, parsed;
	int maxSize = 50;
	int len = 10;
	char* input = (char*)malloc(sizeof(char) * maxSize);
	char* freq = (char*)malloc(sizeof(char) * len);
	int* ptr;

	//deallocation of pre-existing list
	if (freqList != NULL) {
		free(freqList);
		free(freqListSize);
		free(freqIndex);
	}
	
	//new frequency list allocation with default size of 100
	freqListSize = (int*)malloc(sizeof(int));
	*freqListSize = 100;
	freqList = (int*)malloc(sizeof(int)*(*freqListSize));
	freqIndex = (int*)malloc(sizeof(int));
	*freqIndex = 0;

	printf("Frequency Formatting: Single/Multiple separated by ',' Interval connected by '-'\n");
	printf("Frequency: ");
	c = getc(stdin); //flush single char in buffer
	if (isspace(c) ==0) ungetc(c, stdin); //unflush if valid
	while (1){
		c = getc(stdin);
		if (c == EOF || isspace(c)){
			input[i]='\0';
			break;
		}
		input[i]=c;
		if (i == maxSize - 1){
			maxSize *= 2;
			input = (char*)realloc(input, maxSize);
		}
		i++;
	}
	for (j=0; j<=i; j++){
		if (isdigit(input[j])){
			freq[k] = input[j];
			if (k == len -1){
				len *= 2;
				freq = (char*)realloc(freq, len);
			}
			k++;
		} else if (input[j] =='-'){
			freq[k] = '\0'; //insert null terminating char
			parsed = atoi(freq);
			if (parsed > SAMPLE_RATE/2) {
				printf("Error: Frequency %d Hz exceeds the Nyquist Folding Frequency %d\n", parsed, SAMPLE_RATE/2);
				return -1;
			}
			int lowFreq = parsed;
			k = 0;
			for (u = 0; u < (i-j); u++){
				if (input[j+1+u] == ',' || input[j+1+u] == '\0'){
					freq[k] = '\0';
					break;
				} else if (isdigit(input[j+1+u])){
					freq[k] = input[j+1+u];
					k++;
				}
			}
			parsed = atoi(freq);
			if (parsed > SAMPLE_RATE/2) {
				printf("Error: Frequency %d Hz exceeds the Nyquist Folding Frequency %d\n", parsed, SAMPLE_RATE/2);
				return -1;
			} else if (lowFreq >= parsed) {
				printf("Error: Frequency Lower Bound %d Hz equals or exceeds Upper Bound %d Hz\n", lowFreq, parsed);
				return -1;
			}
			int uppFreq = parsed;
			j+=k+1; //advance cursor
			k = 0;
			for (u = lowFreq; u <= uppFreq; u++){
				ptr = (int*)((unsigned char*)freqList + sizeof(int)*(*freqIndex));
				*ptr = u;
				if ((*freqIndex) == (*freqListSize) - 1){
					*freqListSize = (*freqListSize)*2;
					freqList = (int*)realloc(freqList, (*freqListSize)*sizeof(int));
				}
				//printf("Frequency: %d\n", *ptr);			
				*freqIndex = *freqIndex + 1;
			}
		} else if (input[j] ==',' || input[j] == '\0'){
			freq[k] = '\0';
			parsed = atoi(freq);
			ptr = (int*)((unsigned char*)freqList + sizeof(int)*(*freqIndex));
			*ptr = parsed;
			if ((*freqIndex) == (*freqListSize) - 1){
				*freqListSize = (*freqListSize)*2;
				freqList = (int*)realloc(freqList, (*freqListSize)*sizeof(int));
			}
			//printf("Frequency: %d\n", *ptr);
			*freqIndex = *freqIndex + 1;
			k= 0;
		} else {
			continue; //skips over the char
		}
	}

	//sort frequency list via qsort
	qsort(freqList, *(freqIndex), sizeof(int), compare);

	//update all frequency data structure entries
	buildFrequency();
	//deallocate char arrays
	free(input);
	free(freq);

	return 0;
}

//allocation for frequency data structures and thread ids
void buildFrequency(){
	int i;
	destroyFrequency();
	//Float32 compatible
	if (formatType == 0){
		//allocate new memory space for dt
		fq = (Float32DFTSingle**)malloc(sizeof(Float32DFTSingle*) * (*freqIndex));
		for (i=0; i<(*freqIndex); i++){
			fq[i] = (Float32DFTSingle*)malloc(sizeof(Float32DFTSingle));
			fq[i]->realF = (double*)malloc(sizeof(double));
			fq[i]->imagF = (double*)malloc(sizeof(double));
			fq[i]->phaseAngle = (double*)malloc(sizeof(double));
			fq[i]->magnitude = (double*)malloc(sizeof(double));
			fq[i]->decibel = (double*)malloc(sizeof(double));
			fq[i]->zeroIndex = 0;
			fq[i]->frequency = freqList[i];
			fq[i]->samplingRate = SAMPLE_RATE;
		}
	}
	//malloc memory for worker thread ids
	t_id = (pthread_t*)malloc(sizeof(pthread_t)*(*freqIndex));
}

void destroyFrequency(){
	int i;
	//deallocate thread worker ids
	if (t_id !=NULL){
		free(t_id);
	}
	//Float32 compatible
	if (formatType == 0){
		//deallocate previous existing data structure
		if (fq != NULL){
			for (i=0; i<(*freqIndex); i++){
				free(fq[i]->realF);
				free(fq[i]->imagF);
				free(fq[i]->phaseAngle);
				free(fq[i]->magnitude);
				free(fq[i]->decibel);
				free(fq[i]);
			}
			free(fq);
		}
	}
}

void deallocBuffer(){
	//float32 sample type deallocation
	if (formatType == 0){
		if (NUM_CHANNELS == 1){
			free(streamFloat32Mono->singlePhase);
			free(streamFloat32Mono->bufferSize);
			free(streamFloat32Mono->currentTime);
			free(streamFloat32Mono);
		} else if (NUM_CHANNELS ==2) {
			free(streamFloat32Stereo->leftPhase);
			free(streamFloat32Stereo->rightPhase);
			free(streamFloat32Stereo->bufferSize);
			free(streamFloat32Stereo->currentTime);
			free(streamFloat32Stereo);
		}
		free(audioFloat32->inputVolt);
		free(audioFloat32->sampleSize);
		free(audioFloat32->currentTime);
		free(audioFloat32);
	}
}

void DSPdataOutput(){
	int i, format;
	format = __sync_fetch_and_add(&outputFormat, 0); //retrieve value
	if (format == 0){
		//default
		//console version
		for (i=0; i<(*freqIndex); i++){
			printf("Time: "PRINTF_S_FORMAT" Frequency: %d dB: "PRINTF_S_FORMAT"\n", *(audioFloat32->currentTime), fq[i]->frequency, *(fq[i]->decibel));
		}
	} else if (format == 1){
		//gnuPlot version

		/*//temporarily redirect stdout to pipe write-end
		if (dup2(pipeFD[1], STDOUT_FILENO) < 0){
			printf("Error redirecting STDOUT_FILENO to pipe write-end in parent\n");
		} else {
			//redirected region*/

			//signal child of incoming computation packets
			if (pthread_mutex_lock(&(gnuPlotIPC->mutex_IPC)) != 0) {
				//printf("Parent Process Mutex Acquisition Failed\n");
			} else {
				//printf("Parent Process Mutex Acquisition Success\n");

				//update IPC data bundle
				gnuPlotIPC->packetTime = *(audioFloat32->currentTime); //timestamp for each packet
				for (i=0; i<(gnuPlotIPC->packetCount); i++){
					gnuPlotIPC->DSPBundle[i].frequency = fq[i]->frequency;
					gnuPlotIPC->DSPBundle[i].real = *(fq[i]->realF);
					gnuPlotIPC->DSPBundle[i].imag = *(fq[i]->imagF);
					gnuPlotIPC->DSPBundle[i].angle = *(fq[i]->phaseAngle);
					gnuPlotIPC->DSPBundle[i].magnitude = *(fq[i]->magnitude);
					gnuPlotIPC->DSPBundle[i].decibel = *(fq[i]->decibel);
				}
				if (pthread_cond_signal(&(gnuPlotIPC->cond_IPC)) != 0){
					//printf("Parent Process Signal Failed\n");
				} else {
					//printf("Parent Process Signal Success\n");
				}
				if (pthread_mutex_unlock(&(gnuPlotIPC->mutex_IPC)) != 0){
					//printf("Parent Process Mutex Release Failed\n");
				} else {
					//printf("Parent Process Mutex Release Success\n");
				}
			}

			/*//end redirected region
			if (dup2(stdFD[1], STDOUT_FILENO) < 0){
				printf("Error restoring STDOUT_FILENO file descriptor in parent\n");
			}
		}*/
	} else if (format == 2){
		//data file version
		for (i=0; i<(*freqIndex); i++){
			fprintf(fp, PRINTF_S_FORMAT" %d "PRINTF_S_FORMAT"\n", *(audioFloat32->currentTime), fq[i]->frequency, *(fq[i]->decibel));
		}
	} else {

	}
}

void setupPlot(){
	printf("set terminal x11 noraise persist\n");
	printf("unset label\n");
	printf("unset arrow\n");
	printf("set title \"DSP Frequency Decibel RT Spectrum\n");
	printf("set style data lines\n");
	printf("set xlabel \"Time\"\n");
	printf("set ylabel \"Frequency\"\n");
	printf("set zlabel \"Decibel\"\n");
	printf("set xtics 0,1\n");
	printf("set ytics\n");
	printf("set ztics\n");
}

int main(int argc, char const *argv[])
{
	int ret, i, j;
	float* xPtr;
	float* yPtr;
	float* zPtr;
	char ch='x';

	//declare IPC related variables
	pid_t plotterPID; //plotter process PID returned by fork
	FILE* gnuPlotPipe; //gnuPlot stream
	int* graph_Freq; //Y axis
	double* graph_Time; //X axis
	double** graph_Value; //Z axis
	int writePos;

	//allocate shared memory region
	gnuPlotIPC = (IPCSyncObj*)mmap(NULL, sizeof(IPCSyncObj), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
	
	if ((void*)gnuPlotIPC == -1){
		printf("IPC Shared Memory Region Mapping Failed\n");
		return -1;
	}

	gnuPlotIPC->state = 1; //execution switch for plotter 1 = Continue 0 = Exit
	gnuPlotIPC->setup = 0; //setup gnuPlot parameters 0 = uninitialized 1 = initialized
	gnuPlotIPC->interval = 15; //plotting timeframe interval
	gnuPlotIPC->cond_signal = 1; //busy waiting until established buffer size

	if (pthread_condattr_init(&(gnuPlotIPC->condattr_IPC)) != 0 || pthread_condattr_setpshared(&(gnuPlotIPC->condattr_IPC), PTHREAD_PROCESS_SHARED) != 0){
		printf("Error initializing and setting condition variable attribute values\n");
		return -1;
	}

	if (pthread_cond_init(&(gnuPlotIPC->cond_IPC), &(gnuPlotIPC->condattr_IPC)) !=0){
		printf("Error initializing condition variable with attribute\n");
		return -1;
	}

	if (pthread_mutexattr_init(&(gnuPlotIPC->mutexattr_IPC)) != 0 || pthread_mutexattr_setpshared(&(gnuPlotIPC->mutexattr_IPC), PTHREAD_PROCESS_SHARED) != 0){
		printf("Error initializing mutex attribute and set to process shared\n");
		return -1;
	}

	if (pthread_mutex_init(&(gnuPlotIPC->mutex_IPC), &(gnuPlotIPC->mutexattr_IPC)) !=0){
		printf("Error initializing mutex with attributes\n");
		return -1;
	}

	//create pipe in parent before forking
	if (pipe(pipeFD) < 0){
		printf("Failure to create Pipe\n");
		return -1;
	}
	
	//forking
	switch (plotterPID = fork()){
		case -1:
			printf("Failure to fork new child process\n");
			return -1;
		case 0:
			//child process

			printf("Launching...Plotter Started\n");
			//close write end of pipe
			if (close(pipeFD[1]) < 0){
				printf("Failure to close pipe write-end in child\n");
				return -1;
			}
			//store STDOUT_FILENO fd
			stdFD[1] = dup(STDOUT_FILENO);
			if (stdFD[1] < 0){
				printf("Error duplicating STDOUT_FILENO in child\n");
				return -1;
			}
			if (dup2(pipeFD[0], STDIN_FILENO) < 0){
				printf("Error redirecting STDIN_FILENO to pipe read-end in child\n");
				return -1;
			}
			//initialize plotter thread attribute and set joinable state
			if (pthread_attr_init(&child_t_attr) != 0 || pthread_attr_setdetachstate(&child_t_attr, PTHREAD_CREATE_JOINABLE) != 0){
				printf("Error initializing plotter thread attribute and set joinable state in child\n");
				return -1;
			}
			//initialize conditional variable with related attribute and mutex
			if (pthread_condattr_init(&child_condattr_buffer) != 0 || pthread_condattr_setpshared(&child_condattr_buffer, PTHREAD_PROCESS_SHARED) != 0){
				printf("Error initializing and setting condition variable attribute values in child\n");
				return -1;
			}
			if (pthread_cond_init(&child_cond_buffer, &child_condattr_buffer) !=0){
				printf("Error initializing condition variable with attribute in child\n");
				return -1;
			}
			if (pthread_mutexattr_init(&child_mutexattr_buffer) != 0 || pthread_mutexattr_setpshared(&child_mutexattr_buffer, PTHREAD_PROCESS_SHARED) != 0){
				printf("Error initializing mutex attribute and set to process shared in child\n");
				return -1;
			}
			if (pthread_mutex_init(&child_mutex_buffer, &child_mutexattr_buffer) !=0){
				printf("Error initializing mutex with attributes in child\n");
				return -1;
			}

			gnuPlotPipe = popen("gnuplot", "w");
			if (gnuPlotPipe == NULL){
				printf("Failure to open pipe to gnuplot in child\n");
				return -1;
			}
			printf("Launching...gnuPlot Pipe Redirecting\n");
			if (dup2(fileno(gnuPlotPipe), STDOUT_FILENO) < 0){
				printf("Error redirecting STDOUT_FILENO to gnuplot pipe input\n");
				return -1;
			}
			//setup gnuPlot graph appearance if not initialized
			if (__sync_fetch_and_add(&(gnuPlotIPC->setup), 0) != 1){ //retrieve and verify status
				setupPlot(); //proceed to setup
				__sync_fetch_and_add(&(gnuPlotIPC->setup), 1);
			}
			//temporarily busy wait
			while(__sync_fetch_and_add(&(gnuPlotIPC->cond_signal),0)){
				usleep(1); 
			}
			__sync_fetch_and_add(&(gnuPlotIPC->cond_signal),1); //restore value for future use

			//allocate data buffer for gnuPlot
			graph_Freq = (int*)malloc(sizeof(int)*gnuPlotIPC->packetCount);
			graph_Time = (double*)malloc(sizeof(double)*gnuPlotIPC->bufferSize);
			graph_Value = (double**)malloc(sizeof(double*)*gnuPlotIPC->packetCount);
			for (i=0; i<(gnuPlotIPC->packetCount); i++){
				graph_Value[i]=(double*)malloc(sizeof(double)*gnuPlotIPC->bufferSize);
			}

			writePos = 0; //indicate current circular buffer position

			//disable all printf directed to stdout
			while (__sync_fetch_and_add(&(gnuPlotIPC->state), 0) !=0 ){ //atomically verify variable value to continue execution as long as value state != 0
				if (pthread_mutex_lock(&(gnuPlotIPC->mutex_IPC)) != 0) { //acquire lock
					//printf("Child Process Mutex Acquisition Failed\n");
				} else {
					//printf("Child Process Mutex Acquisition Success\n");
					if (pthread_cond_wait(&(gnuPlotIPC->cond_IPC), &(gnuPlotIPC->mutex_IPC)) != 0){
						//printf("Child Process Condition Wait Failed\n");
					} else {
						//printf("Child Process Condition Wait Success\n");
					}
					//write to circular buffer of DSP computation result
					graph_Time[writePos] = gnuPlotIPC->packetTime;						
					for (i=0; i<(gnuPlotIPC->packetCount); i++){
						graph_Freq[i] = gnuPlotIPC->DSPBundle[i].frequency;
						graph_Value[i][writePos] = gnuPlotIPC->DSPBundle[i].decibel;
					}
					writePos = (writePos + 1) % gnuPlotIPC->bufferSize;
					
					if (pthread_mutex_unlock(&(gnuPlotIPC->mutex_IPC)) != 0) { //release lock
						//printf("Child Process Mutex Released Failed\n");
					} else {
						//printf("Child Process Mutex Release Success\n");
					}
				}			
			}
			if (dup2(stdFD[1], STDOUT_FILENO) < 0){
				printf("Error restoring STDOUT_FILENO file descriptor in child\n");
				return -1;
			}
			//destroy plotter thread attribute, mutex attribute, mutex
			if(pthread_attr_destroy(&child_t_attr) != 0 || pthread_mutexattr_destroy(&child_mutexattr_buffer) !=0 || pthread_mutex_destroy(&child_mutex_buffer) !=0){
				printf("Error destroying plotter thread attribute, and/or mutex attribute, and/or mutex in child\n");
				return -1;
			}
			//destroy cond var and associated attributes
			if(pthread_condattr_destroy(&child_condattr_buffer) != 0 || pthread_cond_destroy(&child_cond_buffer) != 0){
				printf("Error destroying conditional variable with associated attribute in child\n");
				return -1;
			}
			printf("Exiting...gnuPlot Realtime Plotting Process Terminating\n");
			free(graph_Freq);
			free(graph_Time);
			for (i=0; i<(gnuPlotIPC->packetCount); i++){
				free(graph_Value[i]);
			}
			free(graph_Value);
			printf("Termination...Child Process Terminated\n");
			//pclose(gnuPlotPipe);
			_Exit(3);
			break;
		default:
			//parent process

			//close read end of pipe
			if (close(pipeFD[0]) < 0){
				printf("Failure to close pipe read-end in parent\n");
				return -1;
			}
			//store STDOUT_FILENO fd
			stdFD[1] = dup(STDOUT_FILENO);
			if (stdFD[1] < 0){
				printf("Error duplicating STDOUT_FILENO in parent\n");
				return -1;
			}

			printf("Launching...Initializing PortAudio\n");
			paErr = Pa_Initialize();
			if (paErr != paNoError) {
				printf("Error Initializing PortAudio\nPortAudio Error: %s",  Pa_GetErrorText(paErr));
				return -1;
			}
			if (queueInputDevice() < 0) return -1;

			//zero both Input/Output Parameters
			memset(&inputParameters, 0, sizeof(PaStreamParameters));
			memset(&outputParameters, 0, sizeof(PaStreamParameters));

			do{		
				ret = selectInput();
			} while (ret < 0);

			do{
				ret = parseFrequency();
			} while (ret < 0);

			//initialize buffer size
			if (*freqIndex > MAX_LOAD) {
				printf("Warning: gnuPlot Plotting Limit Reached, Maximum Handled: %d To Be Generated: %d\n", MAX_LOAD, *freqIndex);
				printf("Warning: gnuPLot Plotting Limit Automatically Adjusted To: %d\n", MAX_LOAD);
				gnuPlotIPC->packetCount = MAX_LOAD;
			} else {
				gnuPlotIPC->packetCount = *freqIndex;
			}
			gnuPlotIPC->bufferSize = gnuPlotIPC->interval * (SAMPLE_RATE/FRAMES_PER_BUFFER);
			__sync_fetch_and_sub(&(gnuPlotIPC->cond_signal),1); //notify child of buffer size calculation completion

			//initialize thread worker attribute and set joinable state
			if (pthread_attr_init(&t_attr) != 0 || pthread_attr_setdetachstate(&t_attr, PTHREAD_CREATE_JOINABLE) != 0){
				printf("Error initializing worker thread attribute and set joinable state\n");
				return -1;
			}

			//initialize cond var attr and set process shared
			if (pthread_condattr_init(&condattr_buffer) != 0 || pthread_condattr_setpshared(&condattr_buffer, PTHREAD_PROCESS_SHARED) != 0){
				printf("Error initializing conditional variable attribute and set to process shared\n");
				return -1;
			}

			//initialize cond var with attribute
			if (pthread_cond_init(&cond_buffer, &condattr_buffer) != 0){
				printf("Error initializing conditional variable\n");
				return -1;
			}

			//initialize mutex attribute and set to process shareable
			if (pthread_mutexattr_init(&mutexattr_buffer) != 0 || pthread_mutexattr_setpshared(&mutexattr_buffer, PTHREAD_PROCESS_SHARED) != 0){
				printf("Error initializing mutex attribute and set to shareable\n");
				return -1;
			}

			//initialize mutex with attribute
			if (pthread_mutex_init(&mutex_buffer, &mutexattr_buffer) != 0){
				printf("Error initializing mutex\n");
				return -1;
			}

			//allocate time variables
			startTime = (struct timespec*)malloc(sizeof(struct timespec));
			streamTime = (struct timespec*)malloc(sizeof(struct timespec));

			//float32 sample type allocation
			if (formatType == 0){
				if (NUM_CHANNELS == 1){	//case mono
					streamFloat32Mono = (activeFloat32Mono*)malloc(sizeof(activeFloat32Mono));
					streamFloat32Mono->singlePhase = (float*)malloc(sizeof(float)*FRAMES_PER_BUFFER);
					streamFloat32Mono->bufferSize = (int*)malloc(sizeof(int));
					streamFloat32Mono->currentTime = (double*)malloc(sizeof(double));
					*(streamFloat32Mono->bufferSize) = FRAMES_PER_BUFFER;
					paErr = Pa_OpenStream(&audioStream, &inputParameters, NULL, SAMPLE_RATE, FRAMES_PER_BUFFER, paNoFlag, streamCallback ,streamFloat32Mono);
				} else if (NUM_CHANNELS == 2){ //case stereo
					streamFloat32Stereo = (activeFloat32Stereo*) malloc(sizeof(activeFloat32Stereo));
					streamFloat32Stereo->leftPhase = (float*)malloc(sizeof(float)*FRAMES_PER_BUFFER);
					streamFloat32Stereo->rightPhase = (float*)malloc(sizeof(float)*FRAMES_PER_BUFFER);
					streamFloat32Stereo->bufferSize = (int*)malloc(sizeof(int));
					streamFloat32Stereo->currentTime = (double*)malloc(sizeof(double));
					*(streamFloat32Stereo->bufferSize) = FRAMES_PER_BUFFER;
					paErr = Pa_OpenStream(&audioStream, &inputParameters, NULL, SAMPLE_RATE, FRAMES_PER_BUFFER, paNoFlag, streamCallback ,streamFloat32Stereo);
				}
			}

			if (paErr != paNoError) {
				printf("Error Opening Stream\nPortAudio Error: %s",  Pa_GetErrorText(paErr));
				return -1;
			}

			//Audio data for continuous processing
			if (formatType == 0){
				//case paFloat32
				audioFloat32 = (DFTFloat32*) malloc(sizeof(DFTFloat32));
				audioFloat32->sampleSize = (int*)malloc(sizeof(int));
				audioFloat32->currentTime = (double*)malloc(sizeof(double));
				if (NUM_CHANNELS == 1){ //case mono
					audioFloat32->inputVolt = (float*)malloc(sizeof(float)*(*streamFloat32Mono->bufferSize));
					*audioFloat32->sampleSize=*streamFloat32Mono->bufferSize;
				} else if (NUM_CHANNELS == 2){ //case stereo
					audioFloat32->inputVolt = (float*)malloc(sizeof(float)*(*streamFloat32Stereo->bufferSize));
					*audioFloat32->sampleSize=*streamFloat32Stereo->bufferSize;
				}
			}
			//printf("*audioFloat32->sampleSize: %d\n", *audioFloat32->sampleSize);
			
			if (outputFormat == 2){//data file output format
				//open data fp
				fp = fopen("Data.txt", "w+");
				if (fp == NULL){
					printf("Error Creating Data File\n");
				}
			}

			paErr = Pa_StartStream(audioStream);
			if (paErr != paNoError) {
				printf("Error Starting Stream\nPortAudio Error: %s",  Pa_GetErrorText(paErr));
				return -1;
			}

			//retrieve stream start time
			if (clock_gettime(CLOCK_MONOTONIC_RAW, startTime) < 0){
				printf("Error Retrieving Stream Start Time\n");
				return -1;
			}

			printf("Running...Stream Is Active\n");	

			//initialize keyboard
			init_keyboard();

			while(ch != 'q') {
				if (formatType == 0){
					//case paFloat32 sample type
					xPtr = (float*) audioFloat32->inputVolt;
					if (NUM_CHANNELS == 1){ //case mono
						yPtr = (float*) streamFloat32Mono->singlePhase;
					} else if (NUM_CHANNELS == 2){ //case stereo
						yPtr = (float*) streamFloat32Stereo->leftPhase;
						zPtr = (float*) streamFloat32Stereo->rightPhase;
					}
					
					//fill missing entries in freq data structure for DFT
					for (i=0; i<(*freqIndex); i++){
						fq[i]->inputSignal = (float*) audioFloat32->inputVolt;
						fq[i]->length = *audioFloat32->sampleSize;
					}

					//begin critical section
					if (pthread_mutex_lock(&mutex_buffer) != 0){ //attempt to lock
						//printf("Main Thread Mutex Acquisition Failed\n"); //failed
					} else {
						//printf("Main Thread Mutex Acquisition Success\n"); //lock acquired
						if (pthread_cond_wait(&cond_buffer, &mutex_buffer) != 0){
							//printf("Main Thread Condition Wait Failed\n");
						} else {
							//printf("Main Thread Condition Wait Success\n");
						}
						if (NUM_CHANNELS == 1){ //case mono
							//copy audio data without processing
							for (i=0; i<*(streamFloat32Mono->bufferSize); i++){
								*xPtr++ = *yPtr++;
							}
							//sync time data
							*(audioFloat32->currentTime) = *(streamFloat32Mono->currentTime);
						} else if (NUM_CHANNELS == 2){
							//equalize stereo to mono
							for (i=0; i<*(streamFloat32Stereo->bufferSize); i++){
								*xPtr++ = (*yPtr++ + *zPtr++)/2.0;
							}
							//sync time data
							*(audioFloat32->currentTime) = *(streamFloat32Stereo->currentTime);
						}
						//end critical section, unlock mutex
						if (pthread_mutex_unlock(&mutex_buffer) != 0){
							//printf("Main Thread Mutex Release Failed\n");
						} else {
							//printf("Main Thread Mutex Release Success\n");
						}
						//dispatch workers
						for (i=0; i<(*freqIndex); i++){
							pthread_create(&(t_id[i]), &t_attr, DFTFloat32SingleFrequency, fq[i]);
						}
						//sync all workers
						for (i=0; i<(*freqIndex); i++){
							pthread_join(t_id[i], NULL);
						}
					}
					//output spectrum result in selected format
					DSPdataOutput();
				}
				if (kbhit()) {
					do{
						ch=readch();
					} while(kbhit());
				}
			}
			//terminate child process
			printf("Exiting...Terminating gnuPlot Realtime Plotting Process\n");
			__sync_fetch_and_sub(&(gnuPlotIPC->state), 1);
			if (pthread_mutex_lock(&(gnuPlotIPC->mutex_IPC)) != 0) {
				//printf("Parent Process Mutex Acquisition Failed\n");
			} else {
				//printf("Parent Process Mutex Acquisition Success\n");
				if (pthread_cond_signal(&(gnuPlotIPC->cond_IPC)) != 0){
					//printf("Parent Process Signal Failed\n");
				} else {
					//printf("Parent Process Signal Success\n");
				}
				//empty signal
				if (pthread_mutex_unlock(&(gnuPlotIPC->mutex_IPC)) != 0){
					//printf("Parent Process Mutex Release Failed\n");
				} else {
					//printf("Parent Process Mutex Release Success\n");
				}
			}
			break;
	}
	close_keyboard();

	printf("Exiting...Closing Stream...\n");
	paErr = Pa_CloseStream(audioStream);
	if (paErr != paNoError) {
		printf("Error Closing Stream\nPortAudio Error: %s\n", Pa_GetErrorText(paErr));
	}

	printf("Exiting...Terminating PortAudio\n");
	paErr = Pa_Terminate();
	if (paErr != paNoError) {
		printf("Error Terminating PortAudio\nPortAudio Error: %s\n", Pa_GetErrorText(paErr));
	}

	//destroy thread worker attribute, mutex attribute, mutex
	if(pthread_attr_destroy(&t_attr) != 0 || pthread_mutexattr_destroy(&mutexattr_buffer) !=0 || pthread_mutex_destroy(&mutex_buffer) !=0){
		printf("Error destroying thread worker attribute, and/or mutex attribute, and/or mutex\n");
		return -1;
	}

	//destroy cond var and associated attributes
	if(pthread_condattr_destroy(&condattr_buffer) != 0 || pthread_cond_destroy(&cond_buffer) != 0){
		printf("Error destroying conditional variable with associated attribute\n");
		return -1;
	}

	if (pthread_condattr_destroy(&(gnuPlotIPC->condattr_IPC)) != 0 || pthread_cond_destroy(&(gnuPlotIPC->cond_IPC)) != 0){
		printf("Error destroying IPC condition variable and associated attributes\n");
		return -1;
	}

	if (pthread_mutexattr_destroy(&(gnuPlotIPC->mutexattr_IPC)) != 0 || pthread_mutex_destroy(&(gnuPlotIPC->mutex_IPC)) != 0){
		printf("Error destroying IPC mutex and associated attributes\n");
		return -1;
	}

	if (munmap(gnuPlotIPC, sizeof(IPCSyncObj)) == -1){
		printf("Unmap Shared Memory Region Failed\n");
		return -1;
	}

	//close data fp
	if (outputFormat == 2 && fp != NULL){
		fclose(fp);
	}

	free(deviceIndex);
	deallocBuffer();
	destroyFrequency();
	free(freqList);
	free(freqListSize);
	free(freqIndex);
	free(startTime);
	free(streamTime);

	printf("Termination...Main Process Terminated\n");
	return 0;
}

void* DFTFloat32SingleFrequency(void* argument){
	//float* inputSignal, double* realF, double* imagF, double* phaseAngle, double* magnitude, double* decibel, int zeroIndex, int frequency, int length, int samplingRate
	int i;
	Float32DFTSingle* DFTData = (Float32DFTSingle* ) argument;
	double realComponent = 0;
	double imagComponent = 0;
	double angle = 0;
	double magn = 0;
	double dB = 0;
	float* ptr = (float*)(DFTData->inputSignal);
	for (i=0; i<DFTData->length; i++){
		realComponent += (*ptr)*cos(2* M_PI * (double)DFTData->frequency/(double)DFTData->samplingRate * (i-(DFTData->zeroIndex)));
		imagComponent += (*ptr)*sin(2* M_PI * (double)DFTData->frequency/(double)DFTData->samplingRate * (i-(DFTData->zeroIndex)));
		ptr++;
	}
	realComponent /= DFTData->length;
	imagComponent /= DFTData->length;
	magn = sqrt(realComponent*realComponent + imagComponent*imagComponent)/(double)DFTData->length;
	dB = 10*log10(magn);
	if (imagComponent == 0 && realComponent == 0){
		angle = 0;
	} else {
		angle = atan2(imagComponent,realComponent)*180.0 / M_PI;
	}
	//printf("realF: "PRINTF_S_FORMAT" imagF: "PRINTF_S_FORMAT" dB: "PRINTF_S_FORMAT" angle: "PRINTF_S_FORMAT"\n", realComponent, imagComponent, dB, angle);
	*(DFTData->realF) = realComponent;
	*(DFTData->imagF) = imagComponent;
	*(DFTData->phaseAngle) = angle;
	*(DFTData->magnitude) = magn;
	*(DFTData->decibel) = dB;
	pthread_exit(0);
}