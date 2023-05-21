#include <ostream>
namespace Color {
    enum Code {
	FG_BLACK    = 30,
        FG_RED      = 31,
        FG_GREEN    = 32,
        FG_YELLOW   = 33,
        FG_BLUE     = 34,
	FG_MAGENTA  = 35,
	FG_CYAN     = 36,
	FG_WHITE    = 37,
        FG_DEFAULT  = 39,
	FG_XBLACK   = 90,
	FG_XRED     = 91,
        FG_XGREEN   = 92,
        FG_XYELLOW  = 93,
        FG_XBLUE    = 94,
	FG_XMAGENTA = 95,
	FG_XCYAN    = 96,
	FG_XWHITE   = 97,

	BG_BLACK    = 40,
        BG_RED      = 41,
        BG_GREEN    = 42,
        BG_YELLOW   = 43,
        BG_BLUE     = 44,
	BG_MAGENTA  = 45,
	BG_CYAN     = 46,
	BG_WHITE    = 47,
        BG_DEFAULT  = 49,
	BG_XBLACK   = 100,
	BG_XRED     = 101,
        BG_XGREEN   = 102,
        BG_XYELLOW  = 103,
        BG_XBLUE    = 104,
	BG_XMAGENTA = 105,
	BG_XCYAN    = 106,
	BG_XWHITE   = 107,

    };
    std::ostream& operator<<(std::ostream& os, Code code) {
        return os << "\033[" << static_cast<int>(code) << "m";
    }
}

