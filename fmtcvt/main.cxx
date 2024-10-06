#include <functional>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <array>
#include <cmath>
#include <tuple>
#include <map>
using namespace std;


// Internal SAO data type.
enum class increment_flag_t : uint8_t {
    SKIP = 0,
    INCREMENT = 1,
    DECREMENT = 2
};


// Convert a character to a `increment_flag_t`.
static inline increment_flag_t c2if(const char _v) {

    const static map<char, increment_flag_t> m = {
        { ' ', increment_flag_t::SKIP },
        { '+', increment_flag_t::INCREMENT },
        { '-', increment_flag_t::DECREMENT }
    };
    return m.at(_v);
}


// Get a magnitude from a spectral type.
// Will throw `out_of_range` if type is composite or invalid.
static inline double st2mag(const array<char, 3>& _type) {

    const static map<char, function<double(const uint8_t)>> m = {
        { 'O', [](const uint8_t _x) -> double { return 0.090 * _x - 5.80; }},
        { 'B', [](const uint8_t _x) -> double { return 0.160 * _x - 4.00; }},
        { 'A', [](const uint8_t _x) -> double { return 0.165 * _x + 0.65; }},
        { 'F', [](const uint8_t _x) -> double { return 0.070 * _x + 2.70; }},
        { 'G', [](const uint8_t _x) -> double { return 0.070 * _x + 4.40; }},
        { 'K', [](const uint8_t _x) -> double { return 0.140 * _x + 5.90; }},
        { 'M', [](const uint8_t _x) -> double { return 0.350 * _x + 8.80; }},
    };
    return m.at(_type[0])(static_cast<uint8_t>(_type[1] == ' ' ? 5 : _type[1] - '0'));
}


// Evaluate space coordinates by RA, DEC and distance.
static tuple<double, double, double> angle2coord(const double _ra, const double _dec, const double _distance) {

    return {
        _distance * cos(_dec) * cos(_ra),
        _distance * cos(_dec) * sin(_ra),
        _distance * sin(_dec)
    };
}


typedef struct __star_coord_t {
    // Distance from the Earth.
    long double distance;
    // Coordinates (where [0, 0, 0] - is a coordinate of the Earth).
    long double x, y, z;
} star_coord_t;


class sao_data final {
public:

    // SAO catalog number - 1-6 bytes number.
    size_t index;

    // Deleted flag.
    bool deleted;

    // RA hours, minutes and seconds.
    uint8_t RAh, RAm;
    double RAs;

    // Annual proper motion and standard deviation in RA, FK4 system.
    double pmRA, epmRA;
    // Time increment.
    increment_flag_t RA2mf;
    // Seconds portion of RA, original epoch precessed to B1950.
    double RA2s;
    // Standard deviation of RA2.
    uint8_t eRA2s;
    // Epoch of RA2 (RA original epoch).
    double EpRA2;
    // Sign Dec, equinox B1950, Epoch=1950.0.
    bool DE;
    // Degrees, minutes and seconds Dec, equinox B1950, Epoch=1950.0.
    uint8_t DEd, DEm;
    double DEs;

    // Annual proper motion in Dec, FK4 system (10).
    double pmDE;
    // Standard deviation of Dec proper motion.
    uint8_t epmDE;
    // DE archminutes flag increment.
    increment_flag_t DE2mf;

    // Seconds of Declination, original epoch, precessed to B1950.
    double DE2s;
    // Standard deviation of DE2.
    uint8_t eDE2s;
    // Epoch of DE2 (Declinaation original epoch).
    double EpDE2;
    // Standard deviation of position at epoch 1950.0.
    uint16_t ePos;
    // Photographic and visual magnitude.
    double Pmag, Vmag;

    // Spectral type, '+++' for composite spectra.
    array<char, 3> SpType;
    // Coded source of visual magnitude; star number and footnotes.
    uint8_t rVmag, rNum;
    // Coded source of photographic magnitude; proper motions; spectral type.
    uint8_t rPmag, rpmRA, rSpType;

    // Coded remarks duplicity and variability.
    uint8_t Rem;

    // Accuracy of V and Ptg: 0 = 2 decimals, 1=1 decimal.
    uint8_t aVmag, aPmag;

    // Source catalog code and number.
    uint8_t rCat;
    uint32_t CatNum;

    // Durchmusterung identification.
    array<char, 13> DM;
    // Henry Draper Catalog (HD or HDE) number.
    array<char, 6> HD;
    // HD component and multiple code.
    uint8_t mHD;
    // Boss General Catalog (GC) number.
    array<char, 5> GC;

    // Right ascension, declination.
    long double RArad, DErad;

    // Hours, minutes and seconds RA, equinox, epoch J2000.0.
    uint8_t RA2000h, RA2000m;
    double RA2000s;
    // Annual proper motion in FK5 system.
    double pmRA2000;
    // Sign Dec, equinox, epoch J2000.0.
    bool DE2000;
    // Degrees, minutes and seconds Dec, equinox, epoch J2000.0.
    uint8_t DE2000d, DE2000m;
    double DE2000s;
    // Annual proper motion in FK5 system (10).
    double pmDE2000;

    // Right ascension and declination, J2000.0, in radians.
    long double RA2000rad, DE2000rad;


    sao_data(const string& _source) {

        this->index = stoull(_source.substr(0, 6));
        this->deleted = _source.at(6) == 'D';
        this->RAh = stoi(_source.substr(7, 2));
        this->RAm = stoi(_source.substr(9, 2));
        this->RAs = stod(_source.substr(11, 6));
        this->pmRA = stod(_source.substr(17, 7));
        this->epmRA = stoi(_source.substr(24, 2));

        this->RA2mf = c2if(_source.at(26));
        this->RA2s = stod(_source.substr(27, 6));
        this->eRA2s = stoi(_source.substr(33, 2));
        this->EpRA2 = stod(_source.substr(35, 6));
        this->DE = _source.at(41) == '+';
        this->DEd = stoi(_source.substr(42, 2));
        this->DEm = stoi(_source.substr(44, 2));
        this->DEs = stod(_source.substr(46, 5));
        this->pmDE = stod(_source.substr(51, 6));
        this->epmDE = stoi(_source.substr(57, 2));

        this->DE2mf = c2if(_source.at(59));
        this->DE2s = stod(_source.substr(60, 5));
        this->eDE2s = stoi(_source.substr(65, 2));
        EpDE2 = stod(_source.substr(67, 6));
        this->ePos = stoul(_source.substr(73, 3));

        this->Pmag = stod(_source.substr(76, 4));
        this->Vmag = stod(_source.substr(80, 4));
        memcpy(&this->SpType[0], &_source.at(84), 3);
        this->rVmag = stoi(_source.substr(87, 2));
        this->rNum = stoi(_source.substr(89, 2));

        this->rPmag = stoi(_source.substr(91, 1));
        this->rpmRA = stoi(_source.substr(92, 1));
        this->rSpType = stoi(_source.substr(93, 1));
        this->Rem = stoi(_source.substr(94, 1));
        this->aVmag = stoi(_source.substr(95, 1));
        this->aPmag = stoi(_source.substr(96, 1));

        this->rCat = stoi(_source.substr(97, 2));
        this->CatNum = stoull(_source.substr(99, 5));

        memcpy(&this->DM[0], &_source.at(104), 13);
        memcpy(&this->HD[0], &_source.at(117), 6);
        this->mHD = _source.at(123);
        memcpy(&this->GC[0], &_source.at(124), 5);

        this->RArad = stold(_source.substr(129, 10));
        this->DErad = stold(_source.substr(139, 11));

        this->RA2000h = stoi(_source.substr(150, 2));
        this->RA2000m = stoi(_source.substr(152, 2));
        this->RA2000s = stod(_source.substr(154, 6));
        this->pmRA2000 = stod(_source.substr(160, 7));

        this->DE2000 = _source.at(167) == '+';
        this->DE2000d = stoi(_source.substr(168, 2));
        this->DE2000m = stoi(_source.substr(170, 2));
        this->DE2000s = stod(_source.substr(172, 5));
        this->pmDE2000 = stod(_source.substr(177, 6));

        this->RA2000rad = stold(_source.substr(183, 10));
        this->DE2000rad = stold(_source.substr(193, 11));
        return;
    }


    // Coordinates evaluated by magnitudes using star's spectral type.
    // Is the star (system) has a composite type, we will not compute its
    // magnitude as for this development step.
    // Throws `invalid_argument` if could not compute the distance.
    star_coord_t get_coordinates(void) const {

        if (this->SpType[0] == ' ' || this->SpType[0] == '+' || this->Pmag == 99.9)
            throw invalid_argument("Could not evaluate distance to the star.");
    
        const double apparent_magniture = 0.5 * this->Pmag + this->Vmag;
        const double absolute_magnitude = st2mag(this->SpType);
        const double distance = pow(10, 0.2 * (apparent_magniture - absolute_magnitude + 5.0));
        const tuple<double, double, double> position = angle2coord(this->RA2000rad, this->DE2000rad, distance);
        
        const star_coord_t coordinates {
            .distance = distance, .x = get<0>(position), .y = get<1>(position), .z = get<2>(position)
        };
        return coordinates;
    }

};


static void convert_sao_to_csv(void) {

    string path;
    cout << "Please prompt a path to a file that contains sao data: ";
    getline(cin, path);

    ofstream csv_file = ofstream("output");
    ifstream sao_file = ifstream(path);
    assert(csv_file.is_open());

    csv_file << "x,y,z,distance,magnitude" << endl;
    if (sao_file.is_open()) {

        string line;
        while (getline(sao_file, line)) try {

            const sao_data data = sao_data(line);
            if (not data.deleted) {

                const star_coord_t coordinates = data.get_coordinates();
                if (coordinates.distance >= 1.0) {

                    csv_file << coordinates.x << ',' << coordinates.y << ',' << coordinates.z << ',';
                    csv_file << coordinates.distance << ',' << st2mag(data.SpType) << endl;
                } else cerr << "[Warning] Skipping star with abnormal distance: " << coordinates.distance << endl;
            } else cout << "[Info] This star had been removed." << endl;
        } catch (const invalid_argument& error) {

            cerr << "[Warning] Could not convert data instance to 3D coordinates. Message:" << endl;
            cerr << " | " << error.what() << endl;
        } catch (const out_of_range& error) {

            cerr << "[Error] An error had been occurred, message:" << endl;
            cerr << " | " << error.what() << endl;
        }
        sao_file.close();
    } else cerr << "[Error] Could not open data file '" << path << "', aborting..." << endl;
    csv_file.close();
    return;
}


int main(void) {

    convert_sao_to_csv();
    return 0;
}
