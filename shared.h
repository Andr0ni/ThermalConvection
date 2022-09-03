//---------------------------------------------------------------------------

#ifndef sharedH
#define sharedH
//---------------------------------------------------------------------------
#include <fstream>
#include <iomanip>
#include <random>
namespace shared
{
    //---------------------------------------------------------------------------
	void cmd(System::UnicodeString s);
    //---------------------------------------------------------------------------
    void cmd(double d);
    //---------------------------------------------------------------------------
    template<typename T>
    void mtrx_write_to_file(
        const T &mtrx, size_t nx, size_t ny, const char* file_name);
    //---------------------------------------------------------------------------
    template<typename T>
    void mtrx_print_to_cnsl(const T mtrx, size_t nx, size_t ny);
    //---------------------------------------------------------------------------
    std::default_random_engine &global_urng()
    {
        static std::default_random_engine u {};
        return u;
    }
    void randomize()
    {
        static std::random_device rd {};
        global_urng().seed(rd());
    }
    int pick(int from, int thru)
    {
        static std::uniform_int_distribution<> d {};
        using parm_t = decltype(d)::param_type;
        return d(global_urng(), parm_t { from, thru });
    }
    double pick(double from, double upto)
    {
        static std::uniform_real_distribution<> d {};
        using parm_t = decltype(d)::param_type;
        return d(global_urng(), parm_t { from, upto });
        //---------------------------------------------------------------------------
    }
} // namespace shared
#endif

