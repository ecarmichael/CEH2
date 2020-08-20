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
#include <Adafruit_NeoPixel.h>

#ifdef __AVR__
#include <avr/power.h>
#endif

// experiment parameters
const int start_stop_pin = 8; // button to start/stop the experiment
int run = 0; //initialize the Run state which is the start of the experiment.


uint32_t ExpTime = 1 * 600000L;

unsigned long tStart = millis();
unsigned long eTime = tStart;

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
//const int Port_R = 7;
// const int Port_C = 8;// from W-maze
//const int Port_L = 9;

//photo beam names
const char* photolabels[] = {"Right", "Right_Center", "Center", "Left_Center", "Left", "Initial"};

// valves
const int solenoidPin_L = 13; // This is the output pin on the Arduino we are using
// int solenoidPin_C = 12; // This is the output pin on the Arduino we are using
const int solenoidPin_R = 11; // This is the output pin on the Arduino we are using
const int valve_L = 0;
// int valve_C = 0;
const int valve_R = 0;
const int valve_dur = 1; // should be in seconds. how long the valve is open.
const int beam_dur = 2; // second before beam break triggers a valve fire.
const int poke_thresh = 3; //how long do they have to nose poke?

// counters
int buttonPushCounter1 = 0;
int buttonPushCounter2 = 0;
int16_t nFlash = 1;
int RTrials = 0;
int LTrials = 0;
int nTrials = 0;
int R_reward = 0;
int L_reward = 0;
int STATE = 6;
int Prior_STATE = 0;
int R_STATE = 0;
int R_first = false;
int L_first = false;

uint16_t i = 0;
int First = 0; // used for printing outputs
int history = 0; // initialize a history of all the states
int prior_state = 0;
//int Loadstate = 0; // toggle state
//int R_beam_state = 0;//toogle state R feeder beam.
int R_oldButtonState = LOW;
int L_oldButtonState = LOW;
int R_toggle = 0;

int R_in = false;
int R_out = true;
int L_in = false;
int L_out = true;

//try a toggle.
int buttonState = 0;
int x = 1;
int first_poke = false;
int beam_time = 0;

// Pixel colors
uint32_t green = pixels.Color(0, 10, 0, 0);
uint32_t red = pixels.Color(10, 0, 0, 0);
uint32_t blue = pixels.Color(0, 0, 10, 0);
uint32_t white = pixels.Color(0, 0, 0, 10);
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

// the loop() methor runs over and over again,
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
    //warp_in();
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

      Serial.print("M_maze");
      Serial.print(';');
      Serial.print("END");
      Serial.print(';');
      Serial.print(buttonPushCounter1);
      Serial.print(';');
      Serial.println(eTime);
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
        tStart = millis();
      }
      else
      {
        digitalWrite(ledPin1, LOW); // set the LED on
        digitalWrite(ledPin2, LOW); // set the LED off
        digitalWrite(ledPinG, LOW); // set the LED off
      }

      // State changes //


      //check that the feeder beams are not broken. Prime the feeders



      //      // Right beam in only
      //      while (digitalRead(beam_R) == LOW)
      //      {
      //        if (first == true)
      //        {
      //         int beam_time = millis();
      //         first = false;
      //         Serial.println("First Break");
      //        }
      //        digitalWrite(ledPin1, HIGH);
      ////        R_in = true;
      ////        R_out = false;
      ////        if (R_in == true)
      ////        {
      ////
      ////          R_in = false; //only use on initial nosepoke.
      //          Prior_STATE = STATE;
      //          //CheckState(beam_R, ledPin1, 0, eTime);
      //          STATE = 0;
      //          delay(10); // prevents gitter in beam breaks.
      ////Serial.print(millis() - beam_time);
      ////Serial.print(';');
      ////Serial.println(R_STATE);
      //          if (((millis() - beam_time) > 3000) && R_STATE == 1)
      //          {
      //                    digitalWrite(ledPin1, LOW);
      //                    delay(100);
      //                            digitalWrite(ledPin1, HIGH);
      //
      //
      //            FireValve(solenoidPin_R, 2);
      //            Serial.print("M_maze");
      //          Serial.print(';');
      //          Serial.print("Right Feeder_fire");
      //          Serial.print(';');
      //          Serial.print(buttonPushCounter1);
      //          Serial.print(';');
      //          Serial.println(eTime);
      //            for (int nFlash = 1; nFlash < 5; nFlash++)
      //            {
      //              // pixels.setPixelColor(0,white, first, count);
      //              pixels.setBrightness(55);
      //              pixels.setPixelColor(STATE, pixels.Color(0, 0, 0, 30));
      //              pixels.show(); // This sends the updated pixel color to the hardware.
      //              delay(100);
      //              pixels.clear();
      //              pixels.show();
      //              delay(100);
      //            }
      //            R_STATE = 0;
      //          }
      //        //}
      //      }



      // RIGHT FEEDER CONTROL
      int R_newButtonState = digitalRead(beam_R);

      // Has the button gone high since we last read it?
      if (R_newButtonState == LOW && R_oldButtonState == HIGH) {
        first_poke = true;
        beam_time = millis();
        buttonPushCounter1++;
        data_log(buttonPushCounter1, "Right Feeder_in", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
        //        Serial.print("M_maze");
        //        Serial.print(';');
        //        Serial.print("Right Feeder_in");
        //        Serial.print(';');
        //        Serial.print(buttonPushCounter1);
        //        Serial.print(';');
        //        Serial.println(eTime);
        digitalWrite(ledPin1, HIGH);

      }

      if (digitalRead(beam_R) == LOW && first_poke == true && R_STATE == true) {
        digitalWrite(ledPin2, LOW);
        //Serial.println(millis() - beam_time);

        if ((millis() - beam_time) > 1000)
        {
          FireValve(solenoidPin_R, 2);
          R_reward++;
          R_STATE = false;
          buttonPushCounter1++;
          data_log(buttonPushCounter1, "Right Feeder_fire", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);

          //          Serial.print("M_maze");
          //          Serial.print(';');
          //          Serial.print("Right Feeder_fire");
          //          Serial.print(';');
          //          Serial.print(buttonPushCounter1);
          //          Serial.print(';');
          //          Serial.println(eTime);
          //          warp_R();
          first_poke = false;
        }

      }
      if (R_newButtonState == HIGH && R_oldButtonState == LOW && beam_time != 0) {
        buttonPushCounter1++;
        data_log(buttonPushCounter1, "Right Feeder_out", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);

        //        Serial.print("M_maze");
        //        Serial.print(';');
        //        Serial.print("Right Feeder_out");
        //        Serial.print(';');
        //        Serial.print(buttonPushCounter1);
        //        Serial.print(';');
        //        Serial.println(eTime);
        first_poke = false;
        digitalWrite(ledPin1, HIGH);

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
        digitalWrite(ledPin2, HIGH);
      }

      if (digitalRead(beam_L) == LOW && first_poke == true && R_STATE == true) {
        digitalWrite(ledPin2, LOW);
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
        digitalWrite(ledPin2, HIGH);
      }
      // Store the button's state so we can tell if it's changed next time round
      L_oldButtonState = L_newButtonState;






      //      if  (digitalRead(beam_R) == HIGH)
      //      {
      //        first = true;
      //      }
      //
      //      if (first == true & digitalRead(beam_R) == LOW)
      //      {
      //        int beam_time = millis();
      //        first = false;
      //        Serial.println("First Break");
      //        while (first == true)
      //        {
      //        if (((millis() - beam_time) > 3000) )
      //        {
      //          FireValve(solenoidPin_R, 2);
      //          Serial.print("M_maze");
      //          Serial.print(';');
      //          Serial.print("Right Feeder_fire");
      //          Serial.print(';');
      //          Serial.print(buttonPushCounter1);
      //          Serial.print(';');
      //          Serial.println(eTime);
      //          for (int nFlash = 1; nFlash < 5; nFlash++)
      //          {
      //            // pixels.setPixelColor(0,white, first, count);
      //            pixels.setBrightness(55);
      //            pixels.setPixelColor(STATE, pixels.Color(0, 0, 0, 30));
      //            pixels.show(); // This sends the updated pixel color to the hardware.
      //            delay(100);
      //            pixels.clear();
      //            pixels.show();
      //            delay(100);
      //          }
      //          R_STATE = 0;
      //        }
      //      }
      //      }

      //      if ((STATE != 0) && (Prior_STATE == 0) && (R_out = false))
      //      {
      //        R_out = true;
      //        R_in = false;
      //        Serial.print("M_maze");
      //        Serial.print(';');
      //        Serial.print("Right Feeder_out");
      //        Serial.print(';');
      //        Serial.print(buttonPushCounter1);
      //        Serial.print(';');
      //        Serial.println(eTime);
      //      }



      //Right feeder beam;  On nosepoke, toggle ON state and start reward counter.

      //      // Get the current state of the button
      //      int newButtonState = digitalRead(beam_R);
      //
      //      // Has the button gone high since we last read it?
      //
      //
      //      // Store the button's state so we can tell if it's changed next time round
      //      oldButtonState = newButtonState;

      //      if (R_beam_state == 0 && digitalRead(beam_R) == HIGH) {
      //        R_beam_state = 1;
      //        Loadstate = !Loadstate;
      //      }
      //      if (R_beam_state == 1 && digitalRead(beam_R) == LOW) {
      //        R_beam_state = 0;
      //      }
      //      if (Loadstate == HIGH) {
      //        // Add Code block
      //        Serial.print("M_maze");
      //        Serial.print(';');
      //        Serial.print("Right Feeder_in");
      //        Serial.print(';');
      //        Serial.print(buttonPushCounter1);
      //        Serial.print(';');
      //        Serial.println(eTime);
      //      }
      //      else {
      //        //Add Code
      //        Serial.print("M_maze");
      //        Serial.print(';');
      //        Serial.print("Right Feeder_out");
      //        Serial.print(';');
      //        Serial.print(buttonPushCounter1);
      //        Serial.print(';');
      //        Serial.println(eTime);
      //      }






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

        pixels.clear();
        pixels.show();
        pixels.setPixelColor(STATE, pixels.Color(0, 10, 0, 0)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
      }
      // Right-center beam only
      while ((digitalRead(beam_CR) == LOW) && STATE != 1)
      {
        Prior_STATE = STATE;
        //CheckState(beam_CR, ledPin1, 1, eTime);
        data_log(buttonPushCounter1, "Right center_beam", nTrials, RTrials, R_reward, LTrials, L_reward, eTime);
        STATE = 1;
      }

      // Left_center beam only
      while ((digitalRead(beam_CL) == LOW) && STATE != 3)
      {
        Prior_STATE = STATE;
        //CheckState(beam_CL, ledPin1, 3, eTime);
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
        pixels.setPixelColor(STATE, pixels.Color(0, 10, 0, 0)); // Moderately bright green color.
        pixels.show(); // This sends the updated pixel color to the hardware.
      }

      //      // Left beam only
      //      while (digitalRead(beam_L) == LOW)
      //      {
      //        Prior_STATE = STATE;
      //        int beam_time = millis();
      //        CheckState(beam_L, ledPin1, 4, eTime);
      //        STATE = 4;
      //
      //        if (((millis() - beam_time) > 3) && R_STATE == true)
      //        {
      //          FireValve(solenoidPin_R, 2);
      //          for (int nFlash = 1; nFlash < 5; nFlash++)
      //          {
      //            // pixels.setPixelColor(0,white, first, count);
      //            pixels.setBrightness(55);
      //            pixels.setPixelColor(STATE, pixels.Color(0, 0, 0, 30));
      //            pixels.show(); // This sends the updated pixel color to the hardware.
      //            delay(100);
      //            pixels.clear();
      //            pixels.show();
      //            delay(100);
      //          }
      //          R_STATE = false;
      //        }
      //      }

      if (digitalRead(Reset_button) == HIGH)
      {

      }
      //exit(0);
    }
  }
  delay(1); // helps prevent false positives.
}

// custom functions

void CheckState(int photo_pin, int LED, int state_id, unsigned long eTime)
{

  //Serial.print(photolabels[state_id]);
  //Serial.print(" Beam broken...");
  digitalWrite(LED, HIGH); // set the LED off
  // delay(500);
  // digitalWrite(ledPin2, LOW);    // set the LED off
  pixels.clear();
  pixels.show();
  pixels.setPixelColor(state_id, pixels.Color(10, 0, 0, 0)); // Moderately bright green color.
  pixels.show(); // This sends the updated pixel color to the hardware.
  First = 1;

  STATE = state_id;
  //Serial.print(".");
  if (STATE == state_id)
  {
    buttonPushCounter1++;
    digitalWrite(LED, LOW);
    //    Serial.println();
    //    Serial.print("This STATE: ");
    //    Serial.print(photolabels[state_id]);
    //    Serial.print("#");
    //    Serial.println(buttonPushCounter1);
    //    Serial.print("From: ");
    //    Serial.print(photolabels[prior_state]);
    //    Serial.print(" to  ");
    //    Serial.print(photolabels[state_id]);
    //    Serial.print(" at ");
    //    Serial.print(eTime / 1000);
    //    Serial.print("s");
    //    Serial.println();
    //    Serial.print("M_maze");
    //    Serial.print(';');
    //    Serial.print(photolabels[state_id]);
    //    Serial.print(';');
    //    Serial.print(buttonPushCounter1);
    //    Serial.print(';');
    //    Serial.println(eTime);

    //
    prior_state = STATE;
    history += STATE;
    STATE = 0;
  }
}


//fire a valve
void FireValve(int valve_pin, int duration)
{
  digitalWrite(valve_pin, HIGH); // Switch Solenoid ON
  delay(duration); // Wait 1 Second
  digitalWrite(valve_pin, LOW); // Switch Solenoid OFF
  //Serial.println("Valve Fired");
}

//log a point in the Maze_data.csv
//// Counter output to serial and csv file.  (from C.D. csBehaviour.ino
//void dataReport(int state_id, unsigned long eTime) {
//  Serial.print("M_maze");
//  Serial.print(';');
//  Serial.print(photolabels[state_id]);
//  Serial.print(';');
//  Serial.print(buttonPushCounter1);
//  Serial.print(';');
//  Serial.println(eTime);
//  //Serial.print(',');
//  //Serial.print(knownValues[0]); //state
//  //Serial.print(',');
//  //Serial.print(knownValues[8]);  //load cell
//  //Serial.print(',');
//  //Serial.print(pulseTrain_chanA[7]); // lick sensor
//  //Serial.print(',');
//  //Serial.print(encoderAngle);     //rotary encoder value
//  //Serial.print(',');
//  //Serial.println(scopeState);
//}
//
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
//void date_str()
//{
//  time_t t = now();
//  Serial.print(year(t));
//  Serial.print("_");
//  Serial.print(month(t));
//  Serial.print("_");
//  Serial.println(day(t));
//}


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

// Log the states
void data_log(int event, const char *ID, int nTrial, int R_trial, int R_reward, int L_trial, int L_reward, int eTime)
{
  Serial.print("M_maze");
  Serial.print(';');
  Serial.print(event);
  Serial.print(';');
  Serial.print(ID);
  Serial.print(';');
  Serial.print(nTrial);
  Serial.print(';');
  Serial.print(R_trial);
  Serial.print(';');
  Serial.print(R_reward);
  Serial.print(';');
  Serial.print(L_trial);
  Serial.print(';');
  Serial.print(L_reward);
  Serial.print(';');
  Serial.println(eTime);
}
//        }
