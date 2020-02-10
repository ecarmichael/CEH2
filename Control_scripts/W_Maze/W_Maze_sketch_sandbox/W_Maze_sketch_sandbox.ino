// Code for 'W' Maze in the Williams lab (McGill).
// Initial version by ecarmichael 2020-02-07 based on Ardunio-AutoReward by 
// jesusjballesteros through open behaviour.  I am very grateful for this framework.  

// Code for the infrared beam-triggered double reward system is shown below.
// Shared first on 05.11.2017 by jesusjballesteros in forum.arduino.cc
// Feel free to modify, share and use it. acknowledgements are always welcome.

#include <Keypad.h>
//neopixel component
#include <Adafruit_NeoPixel.h>
#ifdef __AVR__
  #include <avr/power.h>
#endif

//Main INPUT/OUTPUT system
const int RightValve = 10;          // Right Valve ctrl connected to digital pin 10
const int LeftValve = 9;            // Left Valve ctrl connected to digital pin 13
// photobeams
const int RightLED = 12;            // Left LED anode connected to digital pin 11
const int LeftLED = 11;             // Right LED anode connected to digital pin 12
const int irRightPin = 7;           // Right IR connected to digital pin 4
const int irLeftPin = 8;            // Left IR connected to digital pin 2

//Outputs
#define PIN  6
#define NUMPIXELS 8
Adafruit_NeoPixel pixels = Adafruit_NeoPixel(NUMPIXELS, PIN, NEO_GRB + NEO_KHZ800);

//
byte irRightSensor = 0;             // will read the value from the right IR
byte irLeftSensor = 0;              // will read the value from the left IR
byte alternate = 0;                 // variable to alternate during TRAINING and EXPERIMENT status


//Keyboard implementation (for 8 keys using 5 DIG INPUTS. Can be upgraded)

const byte rows = 2;                 // two rows
const byte cols = 3;                 // three columns
byte rowPins[rows] = {2, 6};         // row pinout of the keypad
byte colPins[cols] = {3, 4, 5};      // column pinout of the keypad
char keys[rows][cols] = {
  {'1', '2', '3',},
  {'4', '5', '6',},
};
Keypad keypad = Keypad( makeKeymap(keys), rowPins, colPins, rows, cols );
char key;
byte menu;
char reset;

//Implementing different STATES. Four states for 3 keys + default.
//¿Possibility of implementing HOLD status to double options without adding new keys?

enum state {HABITUATION, TRAINING, EXPERIMENT1, WAITING, FILLandCLEAN} current_state = WAITING;
      //I define short press 1 = HABITUATION
      //         short press 2 = TRAINING
      //         short press 3 = EXPERIMENT1
      //         short press 4 = FILLandCLEAN
      //         hold 1 = not_def
      //         hold 2 = not_def
      //         hold 3 = not_def


//Initialize values
void setup() {

  Serial.begin(9600);
  keypad.setHoldTime(500);               // Default is 1000mS
  keypad.setDebounceTime(250);           // Default is 50mS

  pinMode(irRightPin, INPUT);        // declare right infrared sensor as input
  pinMode(irLeftPin, INPUT);         // declare left infrared sensor as input
  pinMode(RightLED, OUTPUT);         // declare Right LED pin as output
  pinMode(LeftLED, OUTPUT);          // declare Left LED pin as output
  pinMode(RightValve, OUTPUT);     // declare Right Valve pin as output
  pinMode(LeftValve, OUTPUT);      // declare Left Valve pin as output
  pinMode(LED_BUILTIN, OUTPUT);      // using in_build led (DIG13)

  digitalWrite(irRightPin, HIGH);    // turn on the right IR
  digitalWrite(irLeftPin, HIGH);     // turn on the left IR
  digitalWrite(RightLED, LOW);       // turn off the left LED
  digitalWrite(LeftLED, LOW);        // turn off the left LED

  //neopixels
    pixels.begin(); // This initializes the NeoPixel library.

}

void loop() {

  switch (current_state) {    //Happens depending on the STATE
    case WAITING:
    // pixels.setPixelColor(i, pixels.Color(0,150,0)); // Moderately bright green color.
    //pixels.show(); // This sends the updated pixel color to the hardware.

      digitalWrite(RightLED, HIGH);              //Visual cue for WAITING state
      digitalWrite(LeftLED, HIGH);
      key = keypad.getKey();            // check for key press
      if (key != NO_KEY) {
        menu = key - 48;                // converts ASCII key (char) to menu (byte)
        Serial.println(key);
        switch (menu) {
          case 1:                 // enter the function by short pressing "1"
            {
              current_state = HABITUATION;
              for (byte x = 0; x < 1; x++) {           // Visual cue for HABITUATION state
                digitalWrite(LED_BUILTIN, HIGH);
                delay(500);
                digitalWrite(LED_BUILTIN, LOW);
                delay(500);
              }
            }
            break;

          case 2:                  // enter the function by short pressing "2"
            {
              current_state = TRAINING;
              for (byte x = 0; x < 2; x++) {            // Visual cue for TRAINING state
                digitalWrite(LED_BUILTIN, HIGH);
                delay(500);
                digitalWrite(LED_BUILTIN, LOW);
                delay(500);
              }
            }
            break;

          case 3:                  // enter the function by short pressing "3"
            {
              current_state = EXPERIMENT1;
              for (byte x = 0; x < 3; x++) {            // Visual cue for EXPERIMENT1 state
                digitalWrite(LED_BUILTIN, HIGH);
                delay(500);
                digitalWrite(LED_BUILTIN, LOW);
                delay(500);
              }
            }
            break;

           case 4:                // enter function by short pressing "4"
             {
              current_state = FILLandCLEAN;
              for (byte x = 0; x < 4; x++) {            //Visual cue for FILLandCLEAN state
                digitalWrite(LED_BUILTIN, HIGH);
                delay(500);
                digitalWrite(LED_BUILTIN, LOW);
                delay(500);
                  }
              }
             }
            }
           break;

    case HABITUATION:               // During HABITUATION both IR sensors break generate water delivery, always
      digitalWrite(RightLED, LOW);
      digitalWrite(LeftLED, LOW);
      irRightSensor = digitalRead(irRightPin);      // read right IR value
      irLeftSensor = digitalRead(irLeftPin);        // read left IR value

      if (irLeftSensor == HIGH) {                   // if Left IR does not detect any
        digitalWrite(LeftValve, LOW);               //     Left Valve is OFF
      } else {
        LeftAction();                               // if Left IR detects some, pump ON
      }
      if (irRightSensor == HIGH) {                  // if Left IR does not detect any
        digitalWrite(RightValve, LOW);              //     Left Valve is OFF
      } else {
        RightAction();                              // if Left IR detects some, pump ON
      }

      reset = keypad.getKey();                      // Any key to Reset to WAITING state
      if (reset != NO_KEY) {
        softReset();
      }
      break;

    case TRAINING:                  // During TRAINING alternating IR sensor break delivers water (Starting with LEFT)
      digitalWrite(RightLED, LOW);
      digitalWrite(LeftLED, LOW);
      irRightSensor = digitalRead(irRightPin);       // read right IR value
      irLeftSensor = digitalRead(irLeftPin);         // read left IR value

      if (irLeftSensor == HIGH ) {                   // if Left IR does not detect any
        digitalWrite(LeftValve, LOW);                //     Left Valve is OFF
      } else if ((irLeftSensor == LOW) && (alternate == 0)) {
        LeftAction();                                // if Left IR detects some AND is Left turn, pump ON
        alternate = 1;                               // Changes turn to RIGHT
      }

      if (irRightSensor == HIGH) {                   // if Left IR does not detect any
        digitalWrite(RightValve, LOW);               //      Left Valve is OFF
      } else if ((irRightSensor == LOW) && (alternate == 1)) {
        RightAction();                               // if Left IR detects some AND is Right turn, pump ON
        alternate = 0;                               // Change turn to RIGHT
      }

      reset = keypad.getKey();            // Any key to Reset to WAITING state
      if (reset != NO_KEY) {
        softReset();
      }
      break;

    case EXPERIMENT1:                            // NOT YET DEFINED What happen when State is EXPERIMENT 1
      digitalWrite(RightLED, LOW);
      digitalWrite(LeftLED, LOW);
      irRightSensor = digitalRead(irRightPin);   // read right IR value
      irLeftSensor = digitalRead(irLeftPin);     // read left IR value

      //
      // Here goes code to set EXPERIMENT1 conditions //
      //

      reset = keypad.getKey();            // Any key to Reset to WAITING state
      if (reset != NO_KEY) {
        softReset();
      }
      break;

     case FILLandCLEAN:                 // During FILL/CLEAN IR sensors are not read and Valves are OPEN for 5 secs, every 2 sec
      
      digitalWrite(LeftValve, HIGH);         // Left Valve is turned on
      digitalWrite(LeftLED, HIGH);           // and LED blinks
      digitalWrite(RightValve, HIGH);        // Left Valve is turned on
      digitalWrite(RightLED, HIGH);          // and LED blinks
      delay(5000);
      digitalWrite(LeftValve, LOW);         // Left Valve is turned on
      digitalWrite(LeftLED, LOW);           // and LED blinks
      digitalWrite(RightValve, LOW);        // Left Valve is turned on
      digitalWrite(RightLED, LOW);          // and LED blinks
      delay(2000);
      
      reset = keypad.getKey();            // Any key to Reset to WAITING state
      if (reset != NO_KEY) {
        softReset();
      }
      break; 
    }
}

//Special events

void LeftAction() {                   //IR triggered LeftValve
  digitalWrite(LeftValve, HIGH);          // Left Valve is open briefly
  digitalWrite(LeftLED, HIGH);            // and LED blinks
  delay(100);
  digitalWrite(LeftValve, LOW);
  delay(2000);                            //   for 2 seconds
  digitalWrite(LeftLED, LOW);
  delay(5000);                            // Avoids additional trigger for next 5 sec
}

void RightAction() {                  //IR triggered RightValve
  digitalWrite(RightValve, HIGH);        // Right Valve is turned on
  digitalWrite(RightLED, HIGH);          // and LED blinks
  delay(100);
  digitalWrite(RightValve, LOW);
  delay(2000);                           //   for 2 seconds
  digitalWrite(RightLED, LOW);
  delay(5000);
}

void softReset() {                       // Reset function (taken from arduino.cc forum, user: Volkemon)
  asm volatile ("  jmp 0");
}
