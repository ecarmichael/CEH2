// Rotary Encoder Inputs
#define CLK 2
#define DT 3
#define SW 4

float counter = 0;
int currentStateCLK;
int lastStateCLK;
float abs_pos = 0;
String currentDir = "";
unsigned long lastButtonPress = 0;

void setup() {

  // Set encoder pins as inputs
  pinMode(CLK, INPUT);
  pinMode(DT, INPUT);
  pinMode(SW, INPUT_PULLUP);

  // Setup Serial Monitor
  Serial.begin(9600);

  // Read the initial state of CLK
  lastStateCLK = digitalRead(CLK);
}

void loop() {

  // Read the current state of CLK
  currentStateCLK = digitalRead(CLK);

  // If last and current state of CLK are different, then pulse occurred
  // React to only 1 state change to avoid double count
  if (currentStateCLK != lastStateCLK && currentStateCLK == 1) {

    // If the DT state is different than the CLK state then
    // the encoder is rotating CCW so decrement
    if (digitalRead(DT) != currentStateCLK) {
      counter--;
      currentDir = "CCW";
      if (counter == -24){
        counter = 0;
      }
    } else {
      // Encoder is rotating CW so increment
      counter++;
      currentDir = "CW";
      if (counter == 24){
        counter = 0;
      }
    }

    // convert to degrees of wheel position and catch 0 case
    //if (abs(counter) != 0){
    abs_pos = (abs(counter)/24)*360; 
    //} else {
    //  abs_pos = 0; 
    //}

    Serial.print("Direction: ");
    Serial.print(currentDir);
    Serial.print(" | Counter: ");
    Serial.print(abs(counter));
    Serial.print(" | Rot Deg: ");
    Serial.println(abs_pos);
  }

    //Serial.println(currentStateCLK);

  // Remember last CLK state
  lastStateCLK = currentStateCLK;

  // Read the button state
  int btnState = digitalRead(SW);

  //If we detect LOW signal, button is pressed
  if (btnState == LOW) {
    //if 50ms have passed since last LOW pulse, it means that the
    //button has been pressed, released and pressed again
    if (millis() - lastButtonPress > 50) {
      Serial.println("Button pressed!");
    }

    // Remember last button press event
    lastButtonPress = millis();
  }

  // Put in a slight delay to help debounce the reading
  delay(1);
}