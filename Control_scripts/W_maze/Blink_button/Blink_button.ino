/* LED Blink, Teensyduino Tutorial #1
   http://www.pjrc.com/teensy/tutorial.html
 
   This example code is in the public domain.
*/

// Teensy 2.0 has the LED on pin 11
// Teensy++ 2.0 has the LED on pin 6
// Teensy 3.x / Teensy LC have the LED on pin 13
const int ledPin = 12;
const int buttonPin = 11;

// the setup() method runs once, when the sketch starts

void setup() {
  // initialize the digital pin as an output.
  pinMode(ledPin, OUTPUT);
  Serial.begin(38400);
  pinMode(buttonPin, INPUT);
}

// the loop() methor runs over and over again,
// as long as the board has power

void loop() {
    if (digitalRead(buttonPin) == LOW) {
    Serial.println("Button is not pressed...");
      digitalWrite(ledPin, HIGH);    // set the LED off
      delay(1000);
      digitalWrite(ledPin, LOW);    // set the LED off
  } else {
    Serial.println("Button pressed!!!");
      digitalWrite(ledPin, LOW);   // set the LED on

  }
}
