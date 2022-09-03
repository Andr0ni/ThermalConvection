//---------------------------------------------------------------------------

#pragma hdrstop

#include "shared.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
//---------------------------------------------------------------------------
namespace shared
{
	void cmd(System::UnicodeString s)
    {
        static HANDLE handle;
        if (!handle) {
            AllocConsole();
            handle = GetStdHandle(STD_OUTPUT_HANDLE);
        }
        s += "\n";
        String text = String(s.c_str());

        WriteConsole(handle, text.c_str(), text.Length(), 0, 0);
    }
    //---------------------------------------------------------------------------
    void cmd(double d)
    {
        System::UnicodeString s = FloatToStr(d);
        static HANDLE handle;
        if (!handle) {
            AllocConsole();
            handle = GetStdHandle(STD_OUTPUT_HANDLE);
        }
        s += "\n";
        String text = String(s.c_str());

        WriteConsole(handle, text.c_str(), text.Length(), 0, 0);
    }
    //---------------------------------------------------------------------------
    template<typename T>
    void mtrx_write_to_file(
        const T &mtrx, size_t nx, size_t ny, const char* file_name)
    {
        std::ofstream ofs(file_name);

        if (ofs.is_open()) {
            for (size_t i = 0; i < nx; i++) {
                for (size_t j = 0; j < ny; j++) {
                    ofs << std::fixed << std::setprecision(6) << std::setw(8)
                        << mtrx[i][j] << "\t";
                }
                ofs << std::endl;
            }
        } else {
            std::cout << "writing mtrx to file err";
            exit(-1);
        }
        ofs.close();
    }
    //---------------------------------------------------------------------------
    template<typename T>
    void mtrx_print_to_cnsl(const T mtrx, size_t nx, size_t ny)
    {
        for (size_t i = 0; i < nx; i++) {
            for (size_t j = 0; j < ny; j++) {
                std::cout << std::fixed << std::setprecision(3) << std::setw(8)
                          << mtrx[i][j] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    //---------------------------------------------------------------------------
} // namespace shared

