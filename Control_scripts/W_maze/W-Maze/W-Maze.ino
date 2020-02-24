/* W-Maze Control script,
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


// Teensy 3.x / Teensy LC have the LED on pin 13
const int ledPin1 = 6;
const int ledPin2 = 7;

//neopixel setup
#define NEOPIN         5
#define NUMPIXELS      10
Adafruit_NeoPixel pixels = Adafruit_NeoPixel(NUMPIXELS, NEOPIN, NEO_GRBW + NEO_KHZ800);

const int beam_1 = 2;
const int beam_2 = 3;
const int button_1 = 10;


//valves
int solenoidPin = 11;    //This is the output pin on the Arduino we are using
int valve_L = 0;
int valve_R = 0;

//counters
int buttonPushCounter1 = 0;
int buttonPushCounter2 = 0;
int16_t nFlash = 1;
int STATE = 0;
uint16_t i = 0;

//Pixel colors
uint32_t green = pixels.Color(0, 10, 0, 0);
uint32_t red = pixels.Color(10, 0, 0, 0);
uint32_t blue = pixels.Color(0, 0, 10, 0);
uint32_t white = pixels.Color(00, 0, 0, 10);
uint16_t first = 0;
uint16_t count = NUMPIXELS;

// the setup() method runs once, when the sketch starts

void setup() {
  // initialize the digital pin as an output.
  pinMode(ledPin1, OUTPUT);
  pinMode(ledPin2, OUTPUT);

  //valves
  pinMode(solenoidPin, OUTPUT);           //Sets the pin as an output


  Serial.begin(38400);
  pinMode(beam_1, INPUT);
  pinMode(beam_2, INPUT);
  pinMode(button_1, INPUT);

  pixels.begin(); // This initializes the NeoPixel library.
  pixels.clear();
  pixels.show();
}

// the loop() methor runs over and over again,
// as long as the board has power

void loop() {

  //Rest_button
  if (digitalRead(button_1) == HIGH) {
    Serial.println("RESET Button...");
    //nFlash = 1;
    for (int nFlash = 1; nFlash < 4; nFlash++) {
      digitalWrite(ledPin1, HIGH);    // set the LED off
      digitalWrite(ledPin2, HIGH);    // set the LED off
      pixels.fill(white, first, count);
      //      for (i = 0; i <= 7; i++) pixels.setPixelColor(0, pixels.Color(0, 0, 0, 30));
      pixels.show(); // This sends the updated pixel color to the hardware.

      delay(200);
      digitalWrite(ledPin1, LOW);    // set the LED off
      digitalWrite(ledPin2, LOW);    // set the LED off
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
  } else {
    digitalWrite(ledPin2, LOW);   // set the LED on

  }

  // Right beam

  while (digitalRead(beam_1) == LOW) {
    Serial.println("Left Beam broken...");
    digitalWrite(ledPin1, HIGH);    // set the LED off
    //      delay(500);
    //      digitalWrite(ledPin2, LOW);    // set the LED off
    pixels.clear();
    pixels.show();
    pixels.setPixelColor(7, pixels.Color(10, 0, 0, 0)); // Moderately bright green color.
    pixels.show(); // This sends the updated pixel color to the hardware.

    STATE = 1;
  }
  if (STATE == 1) {
    buttonPushCounter1++;
    digitalWrite(ledPin1, LOW);
    Serial.print("STATE: Right #");
    Serial.println(buttonPushCounter1);

    STATE = 0;
  }
  //    if (digitalRead(beam_1) == LOW) {
  //         buttonPushCounter1++;
  //          Serial.println("Right broken...");
  //      digitalWrite(ledPin1, HIGH);    // set the LED off
  //      delay(500);
  //      digitalWrite(ledPin1, LOW);    // set the LED off
  //
  //      Serial.print("STATE: Right #:  ");
  //      Serial.println(buttonPushCounter1);
  //  } else {
  //      digitalWrite(ledPin1, LOW);   // set the LED on
  //
  //  }

  //beam 2
  while (digitalRead(beam_2) == LOW) {
    Serial.println("Left Beam broken...");
    digitalWrite(ledPin2, HIGH);    // set the LED off
    pixels.clear();
    pixels.show();
    pixels.setPixelColor(1, pixels.Color(10, 0, 0, 0)); // Moderately bright green color.
    pixels.show(); // This sends the updated pixel color to the hardware.

    STATE = 2;
  }
  if (STATE == 2) {
    buttonPushCounter2++;
    digitalWrite(ledPin2, LOW);
    Serial.print("STATE: Left #");
    Serial.println(buttonPushCounter2);
    digitalWrite(solenoidPin, HIGH);    //Switch Solenoid ON
    delay(500);                      //Wait 1 Second
    valve_L++;
    Serial.println("...Reward #: ");
    Serial.println(valve_L);
    digitalWrite(solenoidPin, LOW);     //Switch Solenoid OFF
    STATE = 0;
  }
}
