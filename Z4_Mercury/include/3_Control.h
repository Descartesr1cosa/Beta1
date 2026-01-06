class Control
{
public:
    void Setup(Param *par)
    {
        maximum_Nstep = par->GetInt("max_Nstep");
        tolerance = par->GetDou("tolerance");
        maximum_time = par->GetDou("max_Time");

        file_step = par->GetInt("output_step");
        res_step = par->GetInt("output_residual");
    }

    void UpdateSwitches(const DynamicState &st)
    {
        if_outres = (res_step > 0) && (st.step % res_step == 0);
        if_outfile = (file_step > 0) && (st.step % file_step == 0);

        if_stop = false;
        if_stop = if_stop || (st.step >= maximum_Nstep);
        if_stop = if_stop || (st.time >= maximum_time);
        if_stop = if_stop || (st.residual_max < tolerance);
    }

public:
    bool if_stop = false, if_outres = false, if_outfile = false;

private:
    int maximum_Nstep = 0, file_step = 0, res_step = 0;
    double tolerance = 0.0, maximum_time = 0.0;
};