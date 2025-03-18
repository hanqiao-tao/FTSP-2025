#include "FTSP.h"

bool instance_reader_SALBP_otto(string filename)
{
    ifstream fin;
    fin.open(filename);
    if (fin.is_open() == true)
    {
        string line_str;
        std::getline(fin, line_str);
        assert(line_str == "<number of tasks>");
        std::getline(fin, line_str);
        nj = std::stoi(line_str);

        adj_matr.resize(nj + 1);
        for (int i = 0; i < adj_matr.size(); i++)
        {
            adj_matr[i].resize(nj + 1);
        }

        std::getline(fin, line_str);
        std::getline(fin, line_str);
        assert(line_str == "<cycle time>");
        std::getline(fin, line_str);
        due_date = std::stoi(line_str);

        std::getline(fin, line_str);
        std::getline(fin, line_str);
        assert(line_str == "<order strength>");
        std::getline(fin, line_str);

        std::getline(fin, line_str);
        if (nj != 300) std::getline(fin, line_str);
        std::getline(fin, line_str);
        assert(line_str == "<task times>");

        Job_Time.resize(nj + 1);
        Job_Time[0] = 0;
        for (int j = 1; j <= nj; j++)
        {
            string no_j, dur_j;
            getline(fin, line_str);
            std::stringstream ss(line_str);
            std::getline(ss, no_j, ' ');
            int cur_j = std::stoi(no_j);
            assert(cur_j == j);

            std::getline(ss, dur_j, ' ');
            Job_Time[j] = std::stoi(dur_j);
        }

        std::getline(fin, line_str);
        std::getline(fin, line_str);
        assert(line_str == "<precedence relations>");
        while (getline(fin, line_str) && line_str.empty() == false)
        {
            string pre_s, suc_s;
            std::stringstream ss(line_str);
            std::getline(ss, pre_s, ',');
            int pre = std::stoi(pre_s);
            std::getline(ss, suc_s, ',');
            int suc = std::stoi(suc_s);
            assert(pre < suc);
            adj_matr[pre][suc] = true;
        }
        std::getline(fin, line_str);
        assert(line_str == "<end>");
        fin.close();
        return true;
    }
    else
    {
        std::cerr << "Instance reading Error!!!" << endl;
        return false;
    }
}

bool reader_sort(string filename)
{
    config_matr.push_back(vector<bool>(nj + 1, true));
    ifstream fin;
    fin.open(filename);
    if (fin.is_open() == true)
    {
        string line_str;
        while (getline(fin, line_str) && line_str.empty() == false)
        {
            stringstream ss(line_str);
            string s;
            vector<bool> s_i;
            s_i.push_back(1);
            getline(ss, s, '\040');
            while (s.empty() == false)
            {
                int i = stoi(s);
                if (i == 1)
                {
                    s_i.push_back(1);
                }
                else
                {
                    s_i.push_back(0);
                }
                getline(ss, s, '\040');
                getline(ss, s, '\040');
            }
            config_matr.push_back(s_i);
        }
        fin.close();
        return true;
    }
    else
    {
        cerr << "Sort reading Error!!!" << endl;
        return false;
    }
}

void WriteLogFile(const string& szString, const string& filename)
{
    ofstream fout;
    fout.open(filename, std::ios_base::app);
    fout << szString << endl;
    fout.close();
}

void WriteLogFile_byEntry(const string& szString, const string& filename)
{
    ofstream fout;
    fout.open(filename, std::ios_base::app);
    fout << szString << '\t';
    fout.close();
}

void Writeoutlog(string file_name)
{
    WriteLogFile("", file_name);
    if (if_generate)
    {
        WriteLogFile_byEntry(to_string(LB), file_name);
        WriteLogFile_byEntry(to_string(LB_1), file_name);
        WriteLogFile_byEntry(to_string(LB_machine), file_name);
        WriteLogFile_byEntry(to_string(LB_DW), file_name);
        WriteLogFile_byEntry(to_string(best_solution[0][0]), file_name);
        int ini_sol_0 = BDP_solution[0][0] < H_solution[0][0] ? BDP_solution[0][0] : H_solution[0][0];
        WriteLogFile_byEntry(to_string(ini_sol_0), file_name);
        WriteLogFile_byEntry(to_string(optimalBDP), file_name);
        WriteLogFile_byEntry(to_string(t_lb_1), file_name);
        WriteLogFile_byEntry(to_string(t_machine), file_name);
        WriteLogFile_byEntry(to_string(t_dw), file_name);
        WriteLogFile_byEntry(to_string(t_el), file_name);
        WriteLogFile_byEntry(to_string(t_H), file_name);
        WriteLogFile_byEntry(to_string(t_BDP), file_name);
        WriteLogFile_byEntry(to_string(t_ILS), file_name);
        WriteLogFile_byEntry(to_string(t_pricing), file_name);
    }
    else
    {
        WriteLogFile_byEntry("Fail to find a feasible solution!!!", file_name);
    }
}