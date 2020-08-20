/* M-Maze Control script,
   : reward valves
      _ _ _       _ _ _
     /      \    /      \
    /        \  /        \
   / beam_CL  \/          \ beam_CR
  |                        |
  |                        |
  |                        |
  beam_L * valve_R         beam_R *  valve_L


  Based on sections from:
  Teensyduino Tutorial #1 http://www.pjrc.com/teensy/tutorial.html
  This example code is in the public domain.
  NeoPixel simple sketch (c) 2013 Shae Erisson
  released under the GPLv3 license to match the rest of the AdaFruit NeoPixel library
*/
#include <Adafruit_NeoPixel.h> //runs the neopixel
#ifdef __AVR__
#include <avr/power.h>
#endif
#include <SPI.h>  // for data log
#include <SD.h>

// experiment parameters
int run = 0; //initialize the Run state which is the start of the experiment.
//timing
uint32_t ExpTime = 1 * 600000L;
unsigned long tStart = millis();
unsigned long eTime = tStart;

// ###########Define PINS ####################
const int start_stop_pin = 8; // button to start/stop the experiment
const int Reset_button = 2;
//const int ledPin1 = 4; // used for displaying when the animal nosepokes. 
//const int ledPin2 = 5;

// beams
const int beam_R = 14;
const int beam_CR = 15;
// const int beam_C = 16; // from W-maze
const int beam_CL = 17;
const int beam_L = 18;

//photo beam names
const char* photolabels[] = {"Right", "Right_Center", "Center", "Left_Center", "Left", "Initial"};

// valves
const int solenoidPin_L = 9; // This is the output pin on the Arduino we are using
// int solenoidPin_C = 12; // This is the output pin on the Arduino we are using
const int solenoidPin_R = 9; // This is the output pin on the Arduino we are using
const int valve_L = 0;
// int valve_C = 0;
const int valve_R = 0;
const int valve_dur = 1; // should be in seconds. how long the valve is open.
const int beam_dur = 2; // second before beam break triggers a valve fire.
const int poke_thresh = 3; //how long do they have to nose poke?

// #############STATE Transitions########
int buttonPushCounter1 = 0;
//int16_t nFlash = 1; // for pixel flases
int RTrials = 0;
int LTrials = 0;
int nTrials = 0;
int R_reward = 0;
int L_reward = 0;
int STATE = 6;  //initialize starting state 
int Prior_STATE = 0;
int R_STATE = false;  //init reward state as false. 
int R_first = false;
int L_first = false;

int prior_state = 0;
int R_oldButtonState = LOW;
int L_oldButtonState = LOW;

int R_in = false;
int R_out = true;
int L_in = false;
int L_out = true;

// for nosepokes.
int first_poke = false;
int beam_time = 0;

// neopixel setup
#define NEOPIN 3
#define NUMPIXELS 10
Adafruit_NeoPixel pixels = Adafruit_NeoPixel(NUMPIXELS, NEOPIN, NEO_GRBW + NEO_KHZ800);

uint32_t green = pixels.Color(0, 10, 0, 0); //define some pixel colours
uint32_t red = pixels.Color(10, 0, 0, 0);
uint32_t blue = pixels.Color(0, 0, 10, 0);
uint32_t white = pixels.Color(0, 0, 0, 10);
uint16_t count = NUMPIXELS;


// the setup() method runs once, when the sketch starts
void setup()
{
  // initialize the digital pin as an output.
  //pinMode(ledPin1, OUTPUT);
  //pinMode(ledPin2, OUTPUT);

  // valves
  pinMode(solenoidPin_R, OUTPUT); // Sets the pin as an output
  pinMode(solenoidPin_L, OUTPUT); // Sets the pin as an output

  //beams
  pinMode(beam_R, INPUT);
  pinMode(beam_CR, INPUT);
  //pinMode(beam_C, INPUT);
  pinMode(beam_CL, INPUT);
  pinMode(beam_L, INPUT);

  //buttons
  pinMode(Reset_button, INPUT);
  pinMode(start_stop_pin, INPUT);

  // set reward ports
  //pinMode(Port_R, INPUT);
  //pinMode(Port_C, INPUT);
  //pinMode(Port_L, INPUT);
  pixels.begin(); // This initializes the NeoPixel library.
  pixels.clear();
  pixels.show();

  //serial channels
  // Serial.begin(19200) // used for streaming information to the monitor
  //  Log.begin(9600) // stream the maze state data to a .txt
}

// the loop() method runs over and over again,
// as long as the board has power
void loop()
{

  if (digitalRead(start_stop_pin) == HIGH) // run the control script if the start_stop pin is
  {

    //Serial.println("Starting at time: ");
    if (run == 0)
    {
      run = 255;
    }
    else
    {
      run = 0;
    }

  }
  else
  {
  }
  if (run > 0)
  {
    unsigned long tStart = millis();
    unsigned long eTime = millis() - tStart;

    //print the colum headers
    Serial.println("Task;Event#;STATE;Lap#;R_Lap;R_reward;L_lap;ts(ms)");
    warp_in();
    // end of task print
    if ((millis() - tStart) >= ExpTime) {

        data_log(buttonPushCounter1, "END", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
    }


    // run this code while the counter is less than the ExpTime (experiement time)
    while ((millis() - tStart) < ExpTime) {
      unsigned long eTime = millis() - tStart;
      //Serial.println(eTime);

      // Rest_button
      if (digitalRead(Reset_button) == HIGH)
      {
        Serial.println("RESET Button...");
        // nFlash = 1;
        for (int nFlash = 1; nFlash < 4; nFlash++)
        {
          pixels.setBrightness(44);
          //digitalWrite(ledPin1, HIGH); // set the LED off
          //digitalWrite(ledPin2, HIGH); // set the LED off
          pixels.fill(white, 0, count);
          pixels.show(); // This sends the updated pixel color to the hardware.
          delay(200);
          //digitalWrite(ledPin1, LOW); // set the LED off
          //digitalWrite(ledPin2, LOW); // set the LED off
          pixels.clear();
          pixels.show();
          delay(200);
        }
        buttonPushCounter1 = 0;
        tStart = millis();
      }
      else
      {
        //digitalWrite(ledPin1, LOW); // set the LED on
        //digitalWrite(ledPin2, LOW); // set the LED off
      }

      // State changes //


      // RIGHT FEEDER CONTROL
      int R_newButtonState = digitalRead(beam_R);

      // Has the button gone high since we last read it?
      if (R_newButtonState == LOW && R_oldButtonState == HIGH) {
        first_poke = true;
        beam_time = millis();
        buttonPushCounter1++;
        data_log(buttonPushCounter1, "Right Feeder_in", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
        Pixel_state(STATE, 0,0,0,10); 
        //digitalWrite(ledPin1, HIGH);
      }

      if (digitalRead(beam_R) == LOW && first_poke == true && R_STATE == true) {
//        digitalWrite(ledPin2, LOW);
        if ((millis() - beam_time) > 1000)
        {
          FireValve(solenoidPin_R, 2);
          R_reward++;
          R_STATE = false;
          buttonPushCounter1++;
          data_log(buttonPushCounter1, "Right Feeder_fire", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
          Pixel_state(0, 0,0,10,0); 
          first_poke = false;
        }
      }
      
      if (R_newButtonState == HIGH && R_oldButtonState == LOW && beam_time != 0) {
        buttonPushCounter1++;
        data_log(buttonPushCounter1, "Right Feeder_out", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
        first_poke = false;
        //digitalWrite(ledPin1, HIGH);
      }

      // Store the button's state so we can tell if it's changed next time round
      R_oldButtonState = R_newButtonState;

      //LEFT FEEDER CONTROL
      int L_newButtonState = digitalRead(beam_L);

      // Has the button gone high since we last read it?
      if (L_newButtonState == LOW && L_oldButtonState == HIGH) {
        first_poke = true;
        beam_time = millis();
        buttonPushCounter1++;
        data_log(buttonPushCounter1, "Left Feeder_in", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
        //digitalWrite(ledPin2, HIGH);
      }

      if (digitalRead(beam_L) == LOW && first_poke == true && R_STATE == true) {
        //digitalWrite(ledPin2, LOW);
        if ((millis() - beam_time) > 1000)
        {
          FireValve(solenoidPin_L, 2);
          L_reward++;
          R_STATE = false;
          buttonPushCounter1++;
          data_log(buttonPushCounter1, "Left Feeder_fire", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
          first_poke = false;
        }

      }
      if (L_newButtonState == HIGH && L_oldButtonState == LOW && beam_time != 0) {
        buttonPushCounter1++;
        data_log(buttonPushCounter1, "Left Feeder_out", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
        first_poke = false;
        //digitalWrite(ledPin2, HIGH);
      }
      // Store the button's state so we can tell if it's changed next time round
      L_oldButtonState = L_newButtonState;



      // if the center beam was trigger before the L beam, prime the reward
      if ((Prior_STATE == 3) && STATE == 1)
      {
        if ( R_first == false)
        {
          RTrials++; // count as a right trial.
          nTrials++;
        }
        R_first = true;
        R_STATE = true;
      }

      // Right-center beam only
      while ((digitalRead(beam_CR) == LOW) && STATE != 1)
      {
        Prior_STATE = STATE;
        data_log(buttonPushCounter1, "Right center_beam", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
        STATE = 1;
      }

      // Left_center beam only
      while ((digitalRead(beam_CL) == LOW) && STATE != 3)
      {
        Prior_STATE = STATE;
        data_log(buttonPushCounter1, "Left center_beam", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
        STATE = 3;
      }

      // if the center beam was trigger before the L beam, prime the reward
      if (Prior_STATE == 1 && STATE == 3)
      {
        if ( L_first == false)
        {
          LTrials++;
          nTrials++;
        }
        L_first = true;
        R_STATE = true;
        pixels.clear();
        pixels.show();
        pixels.setPixelColor(STATE, pixels.Color(0, 10, 10, 0)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
      }
      if (R_STATE ==true){
        pixels.clear();
        pixels.show();
        pixels.setPixelColor(STATE, pixels.Color(0, 10, 0, 0)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
      }

      if (digitalRead(Reset_button) == HIGH)
      {
        // do nothing. 
      }
      //exit(0);
    }
  }
  else (
  
  delay(1); // helps prevent false positives.
}

// custom functions

// set the neopixel based on state
void Pixel_state(int STATE, int r, int g, int b, int w){
        pixels.clear();
        pixels.show();
        pixels.setPixelColor(STATE, pixels.Color(r, g, b,w)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
}

//fire a valve
void FireValve(int valve_pin, int duration)
{
  digitalWrite(valve_pin, HIGH); // Switch Solenoid ON
  delay(duration); // Wait 1 Second
  digitalWrite(valve_pin, LOW); // Switch Solenoid OFF
  //Serial.println("Valve Fired");
}


void warp_in() {
  for (int nFlash = 1; nFlash < 5; nFlash++)
  {
    pixels.setBrightness(55);
    pixels.setPixelColor(-1 + nFlash, pixels.Color(0, 0, 0, 30));
    pixels.setPixelColor(8 - nFlash, pixels.Color(0, 0, 0, 30));
    pixels.show(); // This sends the updated pixel color to the hardware.
    delay(100);
    pixels.clear();
    pixels.show();
    delay(100);
  }
}


void warp_out() {
  for (int nFlash = 1; nFlash < 5; nFlash++)
  {
    pixels.setBrightness(55);
    pixels.setPixelColor(3 + nFlash, pixels.Color(0, 0, 0, 30));
    pixels.setPixelColor(4 - nFlash, pixels.Color(0, 0, 0, 30));
    pixels.show(); // This sends the updated pixel color to the hardware.
    delay(100);
    pixels.clear();
    pixels.show();
    delay(100);
  }
}

void warp_R() {
  for (int nFlash = 1; nFlash < 8; nFlash++)
  {
    pixels.setBrightness(55);
    pixels.setPixelColor(7 - nFlash, pixels.Color(0, 0, 0, 30));
    pixels.show(); // This sends the updated pixel color to the hardware.
    delay(20);
    pixels.clear();
    pixels.show();
    delay(20);
  }
}

void warp_L() {
  for (int nFlash = 1; nFlash < 8; nFlash++)
  {
    pixels.setBrightness(55);
    pixels.setPixelColor(0 + nFlash, pixels.Color(0, 0, 0, 30));
    pixels.show(); // This sends the updated pixel color to the hardware.
    delay(20);
    pixels.clear();
    pixels.show();
    delay(20);
  }
}

// Log the states
void data_log(int event, const char *ID, int nTrial, int R_trial, int R_reward, int L_trial, int L_reward, int eTime)
{
  Serial.print("M_maze");
  Serial.print(',');
  Serial.print(event);
  Serial.print(',');
  Serial.print(ID);
  Serial.print(',');
  Serial.print(nTrial);
  Serial.print(',');
  Serial.print(R_trial);
  Serial.print(',');
  Serial.print(R_reward);
  Serial.print(',');
  Serial.print(L_trial);
  Serial.print(',');
  Serial.print(L_reward);
  Serial.print(',');
  Serial.println(eTime);
}

//working toggle state
//      if (newButtonState == HIGH && oldButtonState == LOW) {
//
//        if (R_toggle == 0) {
//          // Toggle on
//            R_toggle  = 1;
//          } else {
//            // Toggle off
//            R_toggle  = 0;
//          }
