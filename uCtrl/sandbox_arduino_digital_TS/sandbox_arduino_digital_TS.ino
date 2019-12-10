//init

#include <TimeLib.h>


//set up pin ids
#define Clockin_pin  2      // Digitial clock in 
#define Clockout_pin  4      // Digitial clock out 

const int ledPin = 13; // onboard LED


// set timing
int Fs = 1000;
float evalEverySample = 1.0;



//void setup() {
//
//
//  analogWriteFrequency(4, 80);
//  analogWrite(4, 128);
//}





void setup()  {
  // set pins initial state
  pinMode(Clockin_pin, INPUT);
  digitalWrite(Clockin_pin, LOW);
  pinMode(Clockout_pin, OUTPUT);
  digitalWrite(Clockout_pin, LOW);
  pinMode(ledPin, OUTPUT);

  // set the Time library to use Teensy 3.0's RTC to keep time
  setSyncProvider(getTeensy3Time);

  Serial.begin(115200);
  while (!Serial);  // Wait for Arduino Serial Monitor to open
  delay(100);
  if (timeStatus() != timeSet) {
    Serial.println("Unable to sync with the RTC");
  } else {
    Serial.println("RTC has set the system time");
  }
}




void loop() {
  for (int x = 2; x < 100; x = x * 1.5) {
  Serial.println(x);

  if x == 42
  
  if (Serial.available()) {
    time_t t = processSyncMessage();
    if (t != 0) {
      Teensy3Clock.set(t); // set the RTC
      setTime(t);
    }
  }
  digitalClockDisplay();
  digitalWrite(ledPin, HIGH);
  delay(500);
  digitalWrite(ledPin, LOW);
  delay(500);

}
}

void digitalClockDisplay() {
  // digital clock display of the time
  Serial.print(hour());
  printDigits(minute());
  printDigits(second());
  Serial.print(" ");
  Serial.print(day());
  Serial.print(" ");
  Serial.print(month());
  Serial.print(" ");
  Serial.print(year());
  Serial.println();
}

time_t getTeensy3Time()
{
  return Teensy3Clock.get();
}

// setup sync message
/*  code to process time sync messages from the serial port   */
#define TIME_HEADER  "T"   // Header tag for serial time sync message

unsigned long processSyncMessage() {
  unsigned long pctime = 0L;
  const unsigned long DEFAULT_TIME = 1357041600; // Jan 1 2013

  if (Serial.find(TIME_HEADER)) {
    pctime = Serial.parseInt();
    return pctime;
    if ( pctime < DEFAULT_TIME) { // check the value is a valid time (greater than Jan 1 2013)
      pctime = 0L; // return 0 to indicate that the time is not valid
    }
  }
  return pctime;
}

void printDigits(int digits) {
  // utility function for digital clock display: prints preceding colon and leading 0
  Serial.print(":");
  if (digits < 10)
    Serial.print('0');
  Serial.print(digits);
}
