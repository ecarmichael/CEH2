//From Arduino to Processing to Txt or cvs etc.
//import
//import processing.serial.*;
//declare
PrintWriter output;
Serial udSerial;

void setup() {
  udSerial = new Serial(this, Serial.list()[0], 9600);
  output = createWriter ("dust_concentration.txt");
}

  void draw() {
    if (udSerial.available() > 0) {
      String SenVal = udSerial.readString();
      if (SenVal != null) {
        output.println(SenVal);
      }
    }
  }

  void keyPressed(){
    output.flush();
    output.close();
    exit();
  }
