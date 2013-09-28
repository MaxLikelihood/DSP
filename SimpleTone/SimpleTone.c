/*

Copyright (c) 2013 Fan Jin under MIT License

*/
//compile:
//gcc -o SimpleTone SimpleTone.c -lportaudio -lm -w

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "portaudio.h"
#include "kbhit.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SAMPLE_FORMAT paFloat32

typedef struct _stream{
	double sampleRate;
	double frequency;
	double volume;
} streamData;
typedef struct _interval{
	unsigned long pause;
	unsigned long length;
} interval;

double playBackFreq;
double deviceSampleRate;
PaError paErr;
streamData audioData;
PaStream* audioStream;
interval audioPattern;
PaStreamParameters inputParameters, outputParameters;

int streamCallback(const void *inputBuffer, void *outputBuffer,unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void *userData){
	unsigned long i;
	float amplitude;
	float* outStream = (float*) outputBuffer;
	streamData* data = (streamData*) userData;

	for (i=0; i<framesPerBuffer; i++){
		amplitude = (float) (sin((double)i * M_PI * 2.0 * ((double)data->frequency / (double)data->sampleRate))) * (double)(data->volume);
		*outStream++ = amplitude;
        *outStream++ = amplitude;
	}
	return paContinue;
}
int selectOutput(){
	int channel;
	PaDeviceIndex index;
	PaDeviceInfo* device;

	printf("Choose Output Device: ");
	scanf("%d", &index);
	if (index >= Pa_GetDeviceCount()) {
		printf("Invalid Index\n");
		return -1;
	}
	device = Pa_GetDeviceInfo(index);
	printf("Output Device Selected: %s\n", device->name);
	printf("Number of Channels: ");
	scanf("%d", &channel);
	outputParameters.device = index;
	outputParameters.channelCount = channel;
	outputParameters.sampleFormat = SAMPLE_FORMAT;
	outputParameters.suggestedLatency = device->defaultLowOutputLatency ;
	outputParameters.hostApiSpecificStreamInfo = NULL;
	deviceSampleRate = device->defaultSampleRate;
	printf("..Verifying Format:\n");
	printf("Device Latency: %f\n", device->defaultLowOutputLatency);
	printf("Sample Rate: %f\n", deviceSampleRate);
	paErr = Pa_IsFormatSupported(NULL, &outputParameters, deviceSampleRate);
	if (paErr == paFormatIsSupported){
		printf("..Format Is Supported\n");
	} else {
		printf("..Format Not Supported\n");
		printf("PortAudio Error: %s\n", Pa_GetErrorText(paErr));
		return -1;
	}
	return 0;
}
int displayDevice(){
	int i;
	PaDeviceIndex deviceCount;
	PaDeviceInfo* currentDevice;
	deviceCount = Pa_GetDeviceCount();
	if (deviceCount >= 0) {
		printf("Device Index\tDevice Name\n");
		for (i = 0; i < deviceCount; i++){
			currentDevice = Pa_GetDeviceInfo(i);
			printf("%d\t\t%s\n", i, currentDevice->name);
		}
		currentDevice = Pa_GetDeviceInfo(Pa_GetDefaultInputDevice());
		printf("Default Input Device: %s\n", currentDevice->name);
		currentDevice = Pa_GetDeviceInfo(Pa_GetDefaultOutputDevice());
		printf("Default Output Device: %s\n", currentDevice->name);
	} else {
		printf("No Device Available\n");
		return -1;
	}
	return 0;
}
void update(){
	double frequency, volume;
	printf("Frequency: ");
	scanf("%lf", &frequency);
	printf("Volume(/1.0): ");
	scanf("%lf", &volume);
	audioData.sampleRate = deviceSampleRate;
	audioData.frequency = frequency;
	audioData.volume = volume;
}
void pattern(){
	unsigned long duration, silence;
	printf("Duration(ms): ");
	scanf("%lu", &duration);
	printf("Pause(ms): ");
	scanf("%lu", &silence);
	audioPattern.length = duration;
	audioPattern.pause = silence;
}
int main(int argc, char const *argv[])
{
	int i, ret;
	char ch='x';

	paErr = Pa_Initialize();
	if (paErr != paNoError) {
		printf("Potential Error Initializing PortAudio\nPortAudio Error: %s\n", Pa_GetErrorText(paErr));
		printf("Exiting...\n");
		return -1;
	}
	printf("Welcome to SimpleTone....\n");
	if (displayDevice() < 0) return -1;
	memset(&inputParameters, 0, sizeof(PaStreamParameters));
	memset(&outputParameters, 0, sizeof(PaStreamParameters));
	do{		
		ret = selectOutput();
	} while (ret < 0);
	update();
	pattern();
	printf("Opening Strean...\n");
	paErr = Pa_OpenStream(&audioStream, NULL, &outputParameters, deviceSampleRate, (deviceSampleRate), paNoFlag, streamCallback, &audioData); //FramePerBuffer adjusted to SampleRate instead of paFramesPerBufferUnspecified
	if (paErr != paNoError) {
		printf("Error Opening Stream\nPortAudio Error: %s\n", Pa_GetErrorText(paErr));
	}
	printf("At any time:\n");
	printf("q - quit\nu - update frequency & volume\np - update duration and pause\n");
	init_keyboard();
	while(ch != 'q') {
		paErr = Pa_StartStream(audioStream);
		if (paErr != paNoError) {
			close_keyboard();
			printf("Error Starting Stream\nPortAudio Error: %s\n", Pa_GetErrorText(paErr));
			init_keyboard();
		}	
		Pa_Sleep(audioPattern.length);	
		paErr = Pa_StopStream(audioStream);
		if (paErr != paNoError) {
			close_keyboard();
			printf("Error Stoping Stream\nPortAudio Error: %s\n", Pa_GetErrorText(paErr));
			init_keyboard();
		}
		if (kbhit()) {
			do{
				ch=readch();
			} while(kbhit());
			puts("");
			if (ch == 'u'){
				close_keyboard();
				update();
				init_keyboard();
			} else if (ch == 'p'){
				close_keyboard();
				pattern();
				init_keyboard();
			}
		}
		if (audioPattern.pause != 0) usleep(audioPattern.pause*1000);
    }	
    printf("Exiting...Closing Stream...\n");
	paErr = Pa_CloseStream(audioStream);
	if (paErr != paNoError) {
		printf("Error Closing Stream\nPortAudio Error: %s\n", Pa_GetErrorText(paErr));
	}
	paErr= Pa_Terminate();
	if (paErr != paNoError) {
		printf("Potential Error Terminating PortAudio\nPortAudio Error: %s\n", Pa_GetErrorText(paErr));
	}
	close_keyboard();
	printf("Terminated\n");
	return 0;
}
