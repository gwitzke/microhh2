/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TIMEDEP
#define TIMEDEP

#include "master.h"
#include "grid.h"
#include "timeloop.h"

class Master;
template<typename> class Grid;

template<typename TF>
class Timedep
{
    public:
        Timedep(Master&, Grid<TF>&, std::string, std::bool);
        ~Timedep();

        bool sw;
        std::string varname;
        std::vector<double> time;
        std::vector<TF> data;


        void create_timedep();
        TF update_time_dependent(Timeloop<TF>&);

        void create_timedep_prof();
        void update_time_dependent_prof(std::vector<TF>, Timeloop<TF>&);

        #ifdef USECUDA
            TF* data_g;
            void update_time_dependent_prof_g(TF*, Timeloop<TF>&);
        #endif

    private:
        Master& master;
        Grid<TF>& grid;
};
#endif