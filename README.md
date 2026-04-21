# Night Sky Creator

Create a beautiful, high-resolution personalized night sky poster for any date, time, and location.

This project generates romantic astronomy-inspired posters with:
- accurate visible star positions
- optional constellation lines
- artistic Milky Way rendering
- nebula and dust lane effects
- high-resolution export for A4 and A3 printing

It is ideal for gifts, anniversaries, birthdays, and special memories.

## Preview

This script creates a poster-style star map with:
- a circular sky chart
- elegant poster framing
- custom title, city, coordinates, and date
- export-ready print files

## Features

- Accurate sky for a chosen **date, time, and location**
- Two rendering styles:
  - `gift` → more artistic and romantic
  - `exact` → more realistic sky appearance
- Optional constellation connection lines
- High-resolution print support:
  - A4
  - A3
- Export formats:
  - PNG
  - PDF
  - TIFF
- Styled Milky Way, nebula, and dust effects
- Clean poster layout for printing or gifting

## Project Structure

text
night_sky_creator/
├─ zvijezde.py
├─ constellationship.fab
├─ de421.bsp
├─ hip_main.dat
├─ README.md
└─ exports/

### File overview

* `zvijezde.py` → main generator script
* `constellationship.fab` → constellation line data
* `de421.bsp` → ephemeris file used for sky calculations
* `hip_main.dat` → star catalog data
* `exports/` → generated poster files

## Installation

Install the required packages:

bash
pip install numpy matplotlib skyfield astropy


## Usage

Run the script:

bash
python zvijezde.py


Generated posters will be saved inside the `exports/` folder.

## Configuration

Edit these values at the top of `zvijezde.py`:

python
MODE="gift"
ENVIRONMENT="dark"

TITLE="START OF SOMETHING NEW"
CITY="TUZLA"
LAT=44.537
LON=18.682

DT_LOCAL=datetime(2026,4,3,22,45,0,tzinfo=ZoneInfo("Europe/Sarajevo"))

PRINT_FORMAT="A3"
PRINT_DPI=600


### Important options

#### Rendering mode

* `gift` → soft artistic look, ideal for presents
* `exact` → closer to realistic sky appearance

#### Environment

* `city`
* `suburban`
* `dark`

#### Print settings

* `A4` or `A3`
* `300 DPI` for testing
* `600 DPI` for high-quality print
* `1200 DPI` for very large and heavy output (can go as high as you want but WILL take more time to render)

## Output

Supported export formats:

* `PNG` → best practical option for printing
* `PDF` → useful when you want vector text and layout
* `TIFF` → optional high-quality print export

Example output filename:

star_poster_gift_A3_600dpi.png
`
## How it works

The script:

1. loads astronomical data
2. computes the visible sky for the chosen location and time
3. projects the sky into a circular poster layout
4. renders stars, constellation lines, Milky Way, nebulae, and dust
5. exports the final print-ready poster

## Best settings for gifts

Recommended setup:

python
MODE="gift"
ENVIRONMENT="dark"
PRINT_FORMAT="A3"
PRINT_DPI=600
EXPORT_PNG=True
EXPORT_PDF=False
EXPORT_TIFF=False

This gives a strong balance between visual quality, file size, and printing reliability.

## Notes

* If you want constellation lines, keep `constellationship.fab` in the same folder as `zvijezde.py`.
* `de421.bsp` is required for accurate sky calculations.
* `hip_main.dat` may be useful for offline use depending on your Skyfield setup.
* For print shops, PNG at 600 DPI is often the most reliable option.

## Example use cases

* anniversary gift
* birthday gift
* first date memory
* engagement present
* astronomy wall art
* custom romantic print

## Future ideas

Possible future improvements:

* GUI version
* live preview mode
* faster rendering
* OpenGL / VisPy renderer
* web app version
* custom color themes
* automatic typography presets

## License

This project is shared for personal and educational use.

If you use or improve it, consider crediting the repository.

## Author

Made with love by Nurka for people who want to turn a special moment into a personalized night sky poster.
