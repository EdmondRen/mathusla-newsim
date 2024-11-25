namespace PARMA
{
    double getHPcpp(int, int, int);
    double getrcpp(double, double);
    double getdcpp(double, double);
    double getSpecCpp(int, double, double, double, double, double);
    double getSpecAngFinalCpp(int, double, double, double, double, double, double);
    double getGenerationCpp(double, double, int *, int, double *, double *);
    double get511fluxCpp(double, double, double);

    class ParmaGen
    {
    public:
        ParmaGen();
    };

}