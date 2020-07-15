/*
 SineWavePoints
 
 Write sine wave points to the serial port, followed by the Carriage Return and LineFeed terminator.
 */

int i = 0;

// The setup routine runs once when you press reset:
void setup() {
  // Initialize serial communication at 9600 bits per second:
  Serial.begin(9600);
}

// The loop routine runs over and over again forever:
void loop() {
  // Write the sinewave points, followed by the terminator "Carriage Return" and "Linefeed".
  Serial.print(sin(i*20.0/360.0));
  Serial.write(13);
  Serial.write(10);
  i += 1;
}
