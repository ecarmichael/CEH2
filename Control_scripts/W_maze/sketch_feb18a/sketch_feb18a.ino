/******************************************************************************
  Template.ino
  A useful starting place when adding a TeensyView to an existing project.

  Marshall Taylor @ SparkFun Electronics, March 15, 2017
  https://github.com/sparkfun/SparkFun_TeensyView_Arduino_Library

  This example sets up the TeensyView and draws a test frame repeatedly.
  The objects in the frame were selected to give copy-paste examples for various
  common operations without a lot of chaff.  See TeensyView.h for specifics.

  Compatible with:
  Teensy LC
  Teensy 3.1
  Teensy 3.2
  Teensy 3.5
  Teensy 3.6

  Development environment specifics:
  Arduino IDE 1.6.12 w/ Teensyduino 1.31
  Arduino IDE 1.8.1 w/ Teensyduino 1.35
  TeensyView v1.0

  This code is released under the [MIT License](http://opensource.org/licenses/MIT).

  Please review the LICENSE.md file included with this example. If you have any questions
  or concerns with licensing, please contact techsupport@sparkfun.com.

  Distributed as-is; no warranty is given.
******************************************************************************/
#include <TeensyView.h>  // Include the SFE_TeensyView library

///////////////////////////////////
// TeensyView Object Declaration //
///////////////////////////////////
//Standard
#define PIN_RESET 15
#define PIN_DC    5
#define PIN_CS    10
#define PIN_SCK   13
#define PIN_MOSI  11

//Alternate (Audio)
//#define PIN_RESET 2
//#define PIN_DC    21
//#define PIN_CS    20
//#define PIN_SCK   14
//#define PIN_MOSI  7


TeensyView oled(PIN_RESET, PIN_DC, PIN_CS, PIN_SCK, PIN_MOSI);

void setup()
{
  oled.begin();    // Initialize the OLED
  oled.clear(ALL); // Clear the display's internal memory
  oled.display();  // Display what's in the buffer (splashscreen)
  delay(1000);     // Delay 1000 ms
  oled.clear(PAGE); // Clear the buffer.

}

void loop()
{
  oled.clear(PAGE);  // Clear the page

  oled.rect(5, 5, 20, 20);  // Draw a rectangle
  oled.rectFill(35, 16, 23, 11);  // Draw a filled rectangle
  oled.circle(22, 20, 7);  // Draw the circle:
  oled.pixel(40, 7, WHITE, NORM);  // Draw a white pixel
  oled.pixel(48, 21, BLACK, NORM);  // Draw a black pixel (on the above rectange)

  oled.setFontType(1);  // Set font to type 1
  oled.setCursor(73, 17); // move cursor
  oled.print("world!");  // Write a byte out as a character
  oled.setFontType(0);  // Set font to type 0
  oled.setCursor(67, 12); // move cursor
  oled.print("Hello");  // Write a byte out as a character

  oled.display();  // Send the PAGE to the OLED memory

  delay(200);
}
