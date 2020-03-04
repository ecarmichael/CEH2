#include <Audio.h>
#include <Wire.h>
#include <SPI.h>
#include <SD.h>
#include <SerialFlash.h>

// GUItool: begin automatically generated code
AudioInputI2S            i2s1;           //xy=241.01041412353516,181.01042556762695
AudioEffectReverb        reverb1;        //xy=504.0104064941406,165.01041412353516
AudioOutputI2S           i2s2;           //xy=765.0104064941406,177.01041412353516
AudioConnection          patchCord1(i2s1, 0, reverb1, 0);
AudioConnection          patchCord2(reverb1, 0, i2s2, 0);
// GUItool: end automatically generated code
