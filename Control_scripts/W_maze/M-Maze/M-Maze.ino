/* M-Maze Control script,
  B_* : beams       R_*: reward valves
      _ _ _       _ _ _
     /      \    /      \
    /    B_CL\  / B_CR   \
  B_L/          \/          \ B_R
  |                        |
  |                        |
  |                        |
  R_left                   R_right


  Based on sections from:
  Teensyduino Tutorial #1 http://www.pjrc.com/teensy/tutorial.html
  This example code is in the public domain.
  NeoPixel simple sketch (c) 2013 Shae Erisson
  released under the GPLv3 license to match the rest of the AdaFruit NeoPixel library
*/
#include <Adafruit_NeoPixel.h>

#ifdef __AVR__
#include <avr/power.h>
#endif

// experiment parameters

uint32_t ExpTime = 1 * 60000L;
uint32_t tStart = millis(); 

const int Reset_button = 2;
// Teensy 3.x / Teensy LC have the LED on pin 13
const int ledPin1 = 4;
const int ledPin2 = 5;
const int ledPinG = 6;

// neopixel setup
#define NEOPIN 3
#define NUMPIXELS 10
Adafruit_NeoPixel pixels = Adafruit_NeoPixel(NUMPIXELS, NEOPIN, NEO_GRBW + NEO_KHZ800);
const int beam_R = 14;
const int beam_CR = 15;
// const int beam_C = 16; // from W-maze
const int beam_CL = 17;
const int beam_L = 18;


// vale sensors
const int Port_R = 7;
// const int Port_C = 8;// from W-maze
const int Port_L = 9;

// valves
int solenoidPin_L = 13; // This is the output pin on the Arduino we are using
// int solenoidPin_C = 12; // This is the output pin on the Arduino we are using
int solenoidPin_R = 11; // This is the output pin on the Arduino we are using
int valve_L = 0;
// int valve_C = 0;
int valve_R = 0;

// counters
int buttonPushCounter1 = 0;
int buttonPushCounter2 = 0;
int16_t nFlash = 1;
int RTrials = 0;
int LTrials = 0;
int nTrials = 0;
int STATE = 0;
int Prior_STATE = 0;
int R_STATE = 0;
uint16_t i = 0;
int First = 0; // used for printing outputs

// Pixel colors
uint32_t green = pixels.Color(0, 10, 0, 0);
uint32_t red = pixels.Color(10, 0, 0, 0);
uint32_t blue = pixels.Color(0, 0, 10, 0);
uint32_t white = pixels.Color(00, 0, 0, 10);
uint16_t first = 0;
uint16_t count = NUMPIXELS;


// the setup() method runs once, when the sketch starts
void setup()
{
  // initialize the digital pin as an output.
  pinMode(ledPin1, OUTPUT);
  pinMode(ledPin2, OUTPUT);
  pinMode(ledPinG, OUTPUT);
  // valves
  pinMode(solenoidPin_R, OUTPUT); // Sets the pin as an output
  pinMode(solenoidPin_L, OUTPUT); // Sets the pin as an output
  Serial.begin(38400);
  pinMode(beam_R, INPUT);
  pinMode(beam_CR, INPUT);
  //pinMode(beam_C, INPUT);
  pinMode(beam_CL, INPUT);
  pinMode(beam_L, INPUT);
  pinMode(Reset_button, INPUT);
  // set reward ports
  pinMode(Port_R, INPUT);
  //pinMode(Port_C, INPUT);
  pinMode(Port_L, INPUT);
  pixels.begin(); // This initializes the NeoPixel library.
  pixels.clear();
  pixels.show();
}

// the loop() methor runs over and over again,
// as long as the board has power
void loop()
{
  while ((millis() - tStart) < ExpTime){
    uint32_t eTime = millis() - tStart; 
    //Serial.print(eTime);   
    // Rest_button
    if (digitalRead(Reset_button) == HIGH)
    {
      Serial.println("RESET Button...");
      // nFlash = 1;
      for (int nFlash = 1; nFlash < 4; nFlash++)
      {
        digitalWrite(ledPin1, HIGH); // set the LED off
        digitalWrite(ledPin2, HIGH); // set the LED off
        digitalWrite(ledPinG, HIGH); // set the LED off
        pixels.fill(white, first, count);
        // for (i = 0; i <= 7; i++) pixels.setPixelColor(0, pixels.Color(0, 0, 0, 30));
        pixels.show(); // This sends the updated pixel color to the hardware.
        delay(200);
        digitalWrite(ledPin1, LOW); // set the LED off
        digitalWrite(ledPin2, LOW); // set the LED off
        pixels.clear();
        pixels.show();
        delay(200);
      }
      buttonPushCounter1 = 0;
      buttonPushCounter2 = 0;
      Serial.print("number of beam 1 breaks:  ");
      Serial.println(buttonPushCounter1);
      Serial.print("number of beam 2 breaks:  ");
      Serial.println(buttonPushCounter2);
    }
    else
    {
      digitalWrite(ledPin1, LOW); // set the LED on
      digitalWrite(ledPin2, LOW); // set the LED off
      digitalWrite(ledPinG, LOW); // set the LED off
    }
    // Right beam
    while (digitalRead(beam_R) == LOW)
    {
      if (First == 0)
      {
        Serial.print("Left Beam broken...");
        digitalWrite(ledPin1, HIGH); // set the LED off
        // delay(500);
        // digitalWrite(ledPin2, LOW);    // set the LED off
        pixels.clear();
        pixels.show();
        pixels.setPixelColor(1, pixels.Color(10, 0, 0, 0)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
        First = 1;
      }
      STATE = 1;
      delay(500);
      Serial.print(".");
    }
    if (STATE == 1)
    {
      buttonPushCounter1++;
      digitalWrite(ledPin1, LOW);
      Serial.println(".");
      Serial.print("STATE: Right #");
      First = 0;
      Prior_STATE = 0;
      STATE = 0;
    }
    // Centre Right beam
    while (digitalRead(beam_CR) == LOW)
    {
      if (First == 0)
      {
        Serial.println();
        Serial.print("Center Right Beam broken...");
        digitalWrite(ledPin1, HIGH); // set the LED off
        digitalWrite(ledPinG, HIGH);
        pixels.clear();
        pixels.show();
        pixels.setPixelColor(2, pixels.Color(10, 0, 0, 0)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
        First = 1;
      }
      STATE = 2;
      delay(500);
      Serial.print(".");
    }
    if (STATE == 2)
    {
      buttonPushCounter1++;
      digitalWrite(ledPin1, LOW);
      digitalWrite(ledPinG, LOW);
      Serial.println();
      Serial.println("STATE: Center Right #");
      First = 0;
      Prior_STATE = 2;
      STATE = 0;
    }
    // Centre beam
    while (digitalRead(beam_R) == LOW)
    {
      if (First == 0)
      {
        Serial.print("Center Beam broken...");
        digitalWrite(ledPinG, HIGH);
        pixels.clear();
        pixels.show();
        pixels.setPixelColor(3, pixels.Color(10, 0, 0, 0)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
        First = 1;
      }
      STATE = 3;
      delay(500);
      Serial.print(".");
    }
    while (digitalRead(Port_R) == HIGH)
    {
      Serial.print("Center Reward...");
      // digitalWrite(ledPinG, HIGH);
      pixels.clear();
      pixels.show();
      pixels.setPixelColor(7, pixels.Color(0, 10, 0, 0)); // Moderately bright green color.
      pixels.show(); // This sends the updated pixel color to the hardware.
      R_STATE = 2;
      delay(200);
      Serial.print(".");
    }
    if (STATE == 3 && R_STATE == 0)
    {
      buttonPushCounter1++;
      digitalWrite(ledPinG, LOW);
      Serial.println();
      Serial.print("STATE: Center #");
      Serial.println(buttonPushCounter1);
      // digitalWrite(solenoidPin_R, HIGH);    //Switch Solenoid ON
      // delay(500);                      //Wait 1 Second
      // valve_L++;
      // Serial.println("...Reward #: ");
      // Serial.println(valve_L);
      // digitalWrite(solenoidPin_R, LOW);     //Switch Solenoid OFF
      nTrials++;
      Serial.print("Correct Right Trials: ");
      Serial.println(RTrials);
      Serial.println("/");
      Serial.println(nTrials);
      Prior_STATE = 3;
      STATE = 0;
    }
    if (Prior_STATE == 3 && R_STATE == 2)
    {
      FireValve(solenoidPin_R, 500);
      for (int nFlash = 1; nFlash < 8; nFlash++)
      {
        pixels.clear();
        pixels.show();
        pixels.setPixelColor(nFlash, pixels.Color(0, 0, 0, 10)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
        delay(50);
        pixels.clear();
        pixels.show();
      }
      pixels.setPixelColor(3, pixels.Color(0, 10, 0, 0)); // Moderately bright green color.
      pixels.show();
      nTrials++;
      Serial.println();
      Serial.print("Correct Right Trials: ");
      Serial.print(RTrials);
      Serial.print("/");
      Serial.println(nTrials);
      Prior_STATE = 3;
      R_STATE = 0;
      STATE = 0;
    }
    else if (R_STATE == 3)
    {
      R_STATE = 0;
    }
    // Centre Left beam
    while (digitalRead(beam_CL) == LOW)
    {
      if (First == 0)
      {
        Serial.print("Center Left Beam broken...");
        digitalWrite(ledPin2, HIGH); // set the LED off
        digitalWrite(ledPinG, HIGH);
        pixels.clear();
        pixels.show();
        pixels.setPixelColor(4, pixels.Color(10, 0, 0, 0)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
        First = 1;
      }
      STATE = 4;
      delay(500);
      Serial.print(".");
    }
    if (STATE == 4)
    {
      buttonPushCounter1++;
      digitalWrite(ledPin2, LOW);
      digitalWrite(ledPinG, LOW);
      Serial.println();
      Serial.print("STATE: Center Left #");
      Serial.println(buttonPushCounter1);
      //
      Prior_STATE = 4;
      STATE = 0;
    }
    // LEFT BEAM
    while (digitalRead(beam_L) == LOW)
    {
      if (First == 0)
      {
        Serial.print("Left Beam broken...");
        digitalWrite(ledPin2, HIGH); // set the LED off
        pixels.clear();
        pixels.show();
        pixels.setPixelColor(5, pixels.Color(10, 0, 0, 0)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
        First = 1;
      }
      STATE = 5;
      delay(500);
      Serial.print(".");
    }
    if (STATE == 5)
    {
      buttonPushCounter2++;
      digitalWrite(ledPin2, LOW);
      Serial.println();
      Serial.print("STATE: Left #");
      Serial.println(buttonPushCounter2);
      // digitalWrite(solenoidPin_L, HIGH);    //Switch Solenoid ON
      // delay(500);                      //Wait 1 Second
      // valve_L++;
      // Serial.println("...Reward #: ");
      // Serial.println(valve_L);
      // digitalWrite(solenoidPin_L, LOW);     //Switch Solenoid OFF
      // Prior_STATE = 5;
      Prior_STATE = 5;
      STATE = 0;
    }
  }
}

// custom functions
void FireValve(int valve_pin, int duration)
{
  digitalWrite(valve_pin, HIGH); // Switch Solenoid ON
  delay(duration); // Wait 1 Second
  digitalWrite(valve_pin, LOW); // Switch Solenoid OFF
  Serial.println("Valve Fired");
}
