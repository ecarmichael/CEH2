#include <Audio.h>
#include <Wire.h>
#include <SPI.h>
#include <SD.h>
#include <SerialFlash.h>

    // Include the audio objects you'll use. For a tone, this is likely a sine wave generator.
#include <AudioSynthSine.h>
#include <AudioOutputI2S.h>

AudioSynthSine     sine1;       // Synthesizer object for the sine wave
AudioOutputI2S     audioOutput; // Output to the headphone/line-out

AudioConnection patchCord1(sine1, 0, audioOutput, 0); // Connect sine wave output to audio output

// pinout

const int t1_pin = 4; // Example: Replace with your desired input pin
const int t2_pin = 6; // Example: Replace with your desired input pin


void setup() {
    // Start serial for debugging messages (optional)
    Serial.begin(9600);
    // Initialize the audio library and start the patch bay
    AudioMemory(8); // Set audio memory size (adjust as needed)
    sine1.frequency(2000); // Set the frequency to A4 (440 Hz)
    sine1.amplitude(0.5); // Set the amplitude (volume) to 0.5 (half volume)
    audioOutput.volume(0.8); // Set the output volume

    // Start the audio object (this will begin playing the tone)
    sine1.begin();
    audioOutput.begin();
}


void loop() {
    // The tone will play continuously because it's not stopped.
    // You can add code to change frequency, amplitude, or stop it here.
    // For example, to play a different note:
    if (digitalRead(t1_pin) == HIGH) {
      sine1.frequency(2000); // Play A5
      sine1.amplitude(0.5); // Set the amplitude (volume) to 0.5 (half volume)
      audioOutput.volume(0.8); // Set the output volume
      sine1.begin();
      audioOutput.begin();
    } else {
      tone1.stop(); // Stop playing if the pin is low
    }

    // Play tone 2 if t2_pin high. 
    if (digitalRead(t2_pin) == HIGH) {
      sine1.frequency(4800); // Play A5
      sine1.amplitude(0.5); // Set the amplitude (volume) to 0.5 (half volume)
      audioOutput.volume(0.8); // Set the output volume
      sine1.begin();
      audioOutput.begin();
    } else {
      tone1.stop(); // Stop playing if the pin is low
    }

    // You can even add a delay to stop the tone for a bit
    // delay(2000);
    // sine1.begin(); // Restart the tone
   }
