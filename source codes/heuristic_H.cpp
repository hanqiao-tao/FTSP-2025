#include "FTSP.h"

vector<int> first_fit_BPP(vector<int> item_list)
{
    int task, remaining_cap;
    int num_item = item_list.size();
    vector<int> ff_sol(num_item, 0);
    if (num_item <= 1) return ff_sol;

    bool over = false;
    int bin_num = 1;
    while (!over)
    {
        remaining_cap = due_date;
        for (int i = 1; i < num_item; i++)
        {
            task = item_list[i];
            if (ff_sol[i] == 0 && Job_Time[task] <= remaining_cap)
            {
                remaining_cap -= Job_Time[task];
                ff_sol[i] = bin_num;
            }
        }

        over = true;
        for (int i = 1; i < num_item; i++)
        {
            if (ff_sol[i] == 0)
            {
                over = false;
                bin_num++;
                break;
            }
        }
    }
    ff_sol[0] = bin_num;
    return ff_sol;
}

void heuristic_H()
{
    clock_t  start_time = clock();
    int L = 0;
    for (int i = 1; i <= nj; i++)
    {
        if (shortest_path[i] > L)
        {
            L = shortest_path[i];
        }
    }
    vector<vector<vector<int>>>task_set_jk(num_config + 1);
    for (int k = 1; k <= num_config; k++)
    {
        task_set_jk[k].resize(L + 1);
        for (int j = 0; j <= L; j++)
        {
            task_set_jk[k][j].reserve(nj + 1);
            task_set_jk[k][j].push_back(0); 
        }
    }
    for (int i = 1; i <= nj; i++)
    {
        int j = shortest_path[i];
        int min_config;
        for (int k = 1; k <= num_config; k++)
        {
            if (config_matr[k][i])
            {
                min_config = k;
                break;
            }
        }
        for (int k = 1; k <= num_config; k++)
        {
            if (config_matr[k][i] && task_set_jk[k][j].size() < task_set_jk[min_config][j].size())
            {
                min_config = k;
            }
        }
        task_set_jk[min_config][j].push_back(i);
    }
    H_solution.resize(2);
    H_solution[0] = vector<int>(nj + 1, 0);
    H_solution[1] = vector<int>(1);
    vector<int> sol_fjk;
    int task_number;
    for (int j = 0; j <= L; j++)
    {
        for (int k = 1; k <= num_config; k++)
        {
            if (task_set_jk[k][j].size() > 1)
            {
                sol_fjk = first_fit_BPP(task_set_jk[k][j]);
                task_number = task_set_jk[k][j].size();
                for (int i = 1; i < task_number; i++)
                {
                    H_solution[0][task_set_jk[k][j][i]] = H_solution[0][0] + sol_fjk[i];
                }
                H_solution[0][0] += sol_fjk[0];
                for (int jj = 0; jj < sol_fjk[0]; jj++)
                {
                    H_solution[1].push_back(k);
                }
            }
        }
    }
    t_H = (double)(clock() - start_time) / CLOCKS_PER_SEC;
}