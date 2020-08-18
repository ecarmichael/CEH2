#include <Audio.h>
#include <Wire.h>
#include <SPI.h>
#include <SD.h>
#include <SerialFlash.h>

// GUItool: begin automatically generated code
AudioInputI2S            i2s1;           //xy=175,128
AudioAnalyzePeak         peak1;          //xy=377,105
AudioAnalyzeRMS          rms1;           //xy=379,154
AudioConnection          patchCord1(i2s1, 0, peak1, 0);
AudioConnection          patchCord2(i2s1, 1, rms1, 0);
// GUItool: end automatically generated code
