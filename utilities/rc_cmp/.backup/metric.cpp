#include "metric.h"


/* metric implemenation --------------------------------------------------*/

metric::metric(points *dataset1, points *dataset2, int var_ds1, int var_ds2,
				ofstream &outfile, string outprefix) : outfile(outfile){

	this->outprefix = outprefix;

	this->dataset1 = dataset1;
	this->dataset2 = dataset2;

	this->var_ds1 = var_ds1;
	this->var_ds2 = var_ds2;
}

void metric::set_dataset1(points *dataset1){
	this->dataset1 = dataset1;
}

void metric::set_dataset2(points *dataset2){
	this->dataset2 = dataset2;
}

string metric::get_res(){
	return string("No Metric Performed");
}

string metric::get_metric_name(){
	return string("Metric");
}

void metric::error_out(string err_ps){

	string out;
	out = "ERROR   |  Metric \"" + this->get_metric_name() + "\":  " + err_ps;
	

	outfile << out  << std::endl;
	cout << out << std::endl;
}



/* mean square error -----------------------------------------------------*/

mse::mse(points *dataset1, points *dataset2, int var_ds1, int var_ds2, 
			ofstream &outfile, string outprefix) 
			: metric(dataset1, dataset2, var_ds1, var_ds2, outfile, outprefix){
	
	return;
}

string mse::get_res(){
	string ret;
	std::stringstream out;
	out.precision(10);
	out << compare();
	return out.str();
}

long double mse::compare(){

	if(!check_card_match())
		return -std::numeric_limits<long double>::infinity();

	long double total_sq_err = 0;

	for(int i=0; i<dataset1->get_num_points(); i++){
		
		pnt p1 = dataset1->get_point(i);
		pnt p2 = dataset2->get_point(i);

        //cerr << var_ds1 << ":" << p1.vals[var_ds1] << " | " << var_ds2 << ":" << p2.vals[var_ds2] << endl;

		if( !isnan(p1.vals[var_ds1]) && !isnan(p2.vals[var_ds2]) ){
			total_sq_err += pow( p1.vals[var_ds1] - p2.vals[var_ds2] , 2 );
		}
	}

	long double res = total_sq_err / ((long double)dataset1->get_num_points());

	return res;
}

bool mse::check_card_match(){	

	if(dataset1->get_num_points() == dataset2->get_num_points())
		return true;

	string err_ps("Mismatched number of coordinates: Dataset 1 = ");
	err_ps += itoa(dataset1->get_num_points());
	err_ps += ", Dataset 2 = ";
	err_ps += itoa(dataset2->get_num_points());

	error_out(err_ps);

	return false;
}

string mse::get_metric_name(){
	return string("Mean Square Error");
}



/* root mean square error ------------------------------------------------*/

rmse::rmse(points *dataset1, points *dataset2, int var_ds1, int var_ds2, 
			ofstream &outfile, string outprefix) 
			: metric(dataset1, dataset2, var_ds1, var_ds2, outfile, outprefix){
	
	return;
}

string rmse::get_res(){
	string ret;
	std::stringstream out;
	out.precision(10);
	out << compare();
	return out.str();
}

long double rmse::compare(){

	if(!check_card_match())
		return -std::numeric_limits<long double>::infinity();

	long double total_sq_err = 0;

	for(int i=0; i<dataset1->get_num_points(); i++){

		pnt p1 = dataset1->get_point(i);
		pnt p2 = dataset2->get_point(i);
	
		if( !isnan(p1.vals[var_ds1]) && !isnan(p2.vals[var_ds2]) )
			total_sq_err += pow( p1.vals[var_ds1] - p2.vals[var_ds2], 2 );
	}

	long double res = sqrt(total_sq_err / ((long double)dataset1->get_num_points()));

	return res;

}

bool rmse::check_card_match(){

	if(dataset1->get_num_points() == dataset2->get_num_points())
		return true;

	string err_ps("Mismatched number of coordinates: Dataset 1 = ");
	err_ps += itoa(dataset1->get_num_points());
	err_ps += ", Dataset 2 = ";
	err_ps += itoa(dataset2->get_num_points());

	error_out(err_ps);

	return false;
}

string rmse::get_metric_name(){
	return string("Root Mean Square Error");
}



/* simple cross correlation / template matching---------------------------*/

scc::scc(points *dataset1, points *dataset2, int var_ds1, int var_ds2,
			ofstream &outfile, string outprefix) 
			: metric(dataset1, dataset2, var_ds1, var_ds2, outfile, outprefix){
	
	return;
}

string scc::get_res(){
	string ret;
	std::stringstream out;
	out.precision(10);
	out << compare();
	return out.str();
}

long double scc::compare(){

	if(!check_card_match())
		return -std::numeric_limits<long double>::infinity();

	long double total = 0;

	for(int i=0; i<dataset1->get_num_points(); i++){

		pnt p1 = dataset1->get_point(i);
		pnt p2 = dataset2->get_point(i);
	
		if( !isnan(p1.vals[var_ds1]) && !isnan(p2.vals[var_ds2]) )
			total += p1.vals[var_ds1] * p2.vals[var_ds2];
	}

	return total;

}

bool scc::check_card_match(){	

	if(dataset1->get_num_points() == dataset2->get_num_points())
		return true;

	string err_ps("Mismatched number of coordinates: Dataset 1 = ");
	err_ps += itoa(dataset1->get_num_points());
	err_ps += ", Dataset 2 = ";
	err_ps += itoa(dataset2->get_num_points());

	error_out(err_ps);

	return false;
}

string scc::get_metric_name(){
	return string("Cross-Correlation / Template Matching");
}



/* map -------------------------------------------------------------------*/

map::map(points *dataset1, points *dataset2, int var_ds1, int var_ds2, 
			ofstream &outfile, string outprefix) 
			: metric(dataset1, dataset2, var_ds1, var_ds2, outfile, outprefix){
	
	return;
}

rgb_color map::get_color_single_sided(long double pt_val, long double range, long double min){
	int val = (int) (((pt_val - min) / range) * ((long double)(255)));

	if(isnan(pt_val)){

		rgb_color nanret;
		nanret.r = 0;
		nanret.g = 255;
		nanret.b = 0;

		return nanret;
	}

	rgb_color ret;
	ret.r = ho_cm[val][0];
	ret.g = ho_cm[val][1];
	ret.b = ho_cm[val][2];

	return ret;

}

rgb_color map::get_color_double_sided(long double pt_val, long double range, long double min,
										long double mid_point){

	if( mid_point >= min && mid_point <= (min+range) ){
	
		long double top_range = (range+min) - mid_point;
		long double bot_range = mid_point - min;

		long double scale_factor;

		if(top_range > bot_range){
			scale_factor = top_range / 127;
		}
		else if(bot_range >= top_range){
			scale_factor = bot_range / 127;
		}	

		int gc_val = (int)(((pt_val-mid_point) / scale_factor) + (long double)(127));

		rgb_color ret;
		ret.r = bty_cm[gc_val][0];
		ret.g = bty_cm[gc_val][1];
		ret.b = bty_cm[gc_val][2];

		return ret;

	}


	return get_color_single_sided(pt_val, range, min);

}



/* difference map --------------------------------------------------------*/

difmap::difmap(indexed_points *dataset1, indexed_points *dataset2, int var_ds1, int var_ds2,
				ofstream &outfile, string outprefix) 
	: map(dataset1, dataset2, var_ds1, var_ds2, outfile, outprefix){	

	return;
}

string difmap::get_res(){
	return compare();
}

string difmap::compare(){

	indexed_points *ds1_ip = dynamic_cast<indexed_points*>(dataset1);
	indexed_points *ds2_ip = dynamic_cast<indexed_points*>(dataset2);

	if(ds1_ip == NULL || ds2_ip == NULL)
		return string("Improper Datasets");

	string ret("");

	if(!check_card_match()){
		ret = "No Output";
		return ret;
	}

	long double max = get_max_dif();
	long double min = get_min_dif();

	long double range = fabs( max - min );

	CImg<unsigned char> imgOut( ds1_ip->get_dim().sizes[0], 
								ds1_ip->get_dim().sizes[1],
								ds1_ip->get_dim().sizes[2], 
								3, 0 );	

	for(int x=0; x<ds1_ip->get_dim().sizes[0]; x++){
		for(int y=0; y<ds1_ip->get_dim().sizes[1]; y++){
			for(int z=0; z<ds1_ip->get_dim().sizes[2]; z++){

				//TODO: loop over every dep_var

				layout loc;
				loc.arr_size = 3;
				loc.sizes = new int[3];
				loc.sizes[0] = x;
				loc.sizes[1] = y;
				loc.sizes[2] = z;

				long double pt_val = fabs( ds1_ip->get_indexed_point(loc).vals[var_ds1] 
										   - ds2_ip->get_indexed_point(loc).vals[var_ds2] );

				rgb_color col = get_color_double_sided(pt_val, range, min, 0);

				imgOut(x, ds1_ip->get_dim().sizes[1] - y - 1, 0, 0) = col.r;
				imgOut(x, ds1_ip->get_dim().sizes[1] - y - 1, 0, 1) = col.g;
				imgOut(x, ds1_ip->get_dim().sizes[1] - y - 1, 0, 2) = col.b;
			}
		}
	}

	//Determine if one dataset is always larger
	int large_ds;
	if( 0 >= min && 0 <= max )
		large_ds = -1;
	else{

		layout loc;
		loc.arr_size = 3;
		loc.sizes = new int[3];
		loc.sizes[0] = 0;
		loc.sizes[1] = 0;
		loc.sizes[2] = 0;

		if(ds1_ip->get_indexed_point(loc).vals[var_ds1] 
			> ds2_ip->get_indexed_point(loc).vals[var_ds2]){

			large_ds = 1;
		}
		else
			large_ds = 2;
	}

	//Construct file name
	string output_name(outprefix);
	output_name += "_difmap_";
	output_name += "var1_";
	output_name += itoa(var_ds1);
	output_name += "_var2_";
	output_name += itoa(var_ds2);

	//Write info to txt file
	output_info(output_name + ".txt", max, min, large_ds);

	//Write image
	imgOut.save_bmp((output_name + ".bmp").c_str());

	ret = output_name;
	return ret;
}

bool difmap::check_card_match(){	
	if(dataset1->get_num_points() != dataset2->get_num_points()){

		string err_ps("Mismatched number of coordinates: Dataset 1 = ");
		err_ps += itoa(dataset1->get_num_points());
		err_ps += ", Dataset 2 = ";
		err_ps += itoa(dataset2->get_num_points());

		error_out(err_ps);
		
		return false;
	}

	//TODO
	//Check dimensions

	return true;
}

string difmap::get_metric_name(){
	return string("Difference Map");
}

long double difmap::get_max_value(int ds){

	long double max = -std::numeric_limits<long double>::infinity();

	int size = dataset1->get_num_points();

	for(int i=0; i<size; i++){

		if(ds == 1 && dataset1->get_point(i).vals[var_ds1] > max){
			max = dataset1->get_point(i).vals[var_ds1];
		}
		if(ds == 2 && dataset2->get_point(i).vals[var_ds2] > max){
			max = dataset2->get_point(i).vals[var_ds2];
		}
	}

	return max;
}

long double difmap::get_min_value(int ds){

	long double min = std::numeric_limits<long double>::infinity();

	int size = dataset1->get_num_points();

	for(int i=0; i<size; i++){

		if(ds == 1 && dataset1->get_point(i).vals[var_ds1] < min){
			min = dataset1->get_point(i).vals[var_ds1];
		}
		if(ds == 2 && dataset2->get_point(i).vals[var_ds2] < min){
			min = dataset2->get_point(i).vals[var_ds2];
		}

	}

	return min;
}


long double difmap::get_max_dif(){

	long double maxdif = -std::numeric_limits<long double>::infinity();

	int size = dataset1->get_num_points();

	for(int i=0; i<size; i++){
		
		long double dif = fabs( dataset1->get_point(i).vals[var_ds1] 
								- dataset2->get_point(i).vals[var_ds2] );
		
		if(dif > maxdif)
			maxdif = dif;
	}

	return maxdif;
}

long double difmap::get_min_dif(){

	long double mindif = std::numeric_limits<long double>::infinity();

	int size = dataset1->get_num_points();

	for(int i=0; i<size; i++){

		long double dif = fabs( dataset1->get_point(i).vals[var_ds1] 
								- dataset2->get_point(i).vals[var_ds2] );
		if(dif < mindif)
			mindif = dif;
	}

	return mindif;
}

void difmap::output_info(string filename, long double max, long double min, int large_ds){
	ofstream infowrite;
	infowrite.open(filename.c_str(), ofstream::out);

    infowrite << "FILE 1 MAX: " << get_max_value(1) << endl;
    infowrite << "FILE 1 MIN: " << get_min_value(1) << endl;

    infowrite << "FILE 2 MAX: " << get_max_value(2) << endl;
    infowrite << "FILE 2 MIN: " << get_min_value(2) << endl;

	infowrite << "Max Difference: " << max << endl;
	infowrite << "Min Difference: " << min << endl;
	infowrite << "Range: " << (max-min) << endl;

	if(large_ds != -1)
		infowrite << "Dataset with largest values: " << large_ds << endl;

	infowrite << "\n" << "Light green points: Values that could not be interpolated" << endl;
	infowrite << "Dark green points: Divider" << endl;
}



/* difmap with key -------------------------------------------------------*/

difmap_wkey::difmap_wkey(indexed_points *dataset1, indexed_points *dataset2, 
							int var_ds1, int var_ds2,
							ofstream &outfile, string outprefix) 
	: difmap(dataset1, dataset2, var_ds1, var_ds2, outfile, outprefix){	

	border_width = 4;
	key_width = 20;

	border_color.r = 0;
	border_color.g = 128;
	border_color.b = 0;

	return;
}

string difmap_wkey::get_metric_name(){
	return string("Difference Map with Key");
}

string difmap_wkey::compare(){

	indexed_points *ds1_ip = dynamic_cast<indexed_points*>(dataset1);
	indexed_points *ds2_ip = dynamic_cast<indexed_points*>(dataset2);

	if(ds1_ip == NULL || ds2_ip == NULL)
		return string("Improper Datasets");

	string ret("");

	if(!check_card_match()){
		ret = "No Output";
		return ret;
	}

	long double max = get_max_dif();
	long double min = get_min_dif();

	long double range = fabs( max - min );

	CImg<unsigned char> imgOut( ds1_ip->get_dim().sizes[0] + border_width + key_width, 
								ds1_ip->get_dim().sizes[1],
								ds1_ip->get_dim().sizes[2], 
								3, 0 );	

	//Write difference map to image
	for(int x=0; x<ds1_ip->get_dim().sizes[0]; x++){
		for(int y=0; y<ds1_ip->get_dim().sizes[1]; y++){
			for(int z=0; z<ds1_ip->get_dim().sizes[2]; z++){
				//TODO: loop over every dep_var

				layout loc;
				loc.arr_size = 3;
				loc.sizes = new int[3];
				loc.sizes[0] = x;
				loc.sizes[1] = y;
				loc.sizes[2] = z;

				long double pt_val = fabs( ds1_ip->get_indexed_point(loc).vals[var_ds1]
										   - ds2_ip->get_indexed_point(loc).vals[var_ds2] );

				rgb_color col = get_color_double_sided(pt_val, range, min, 0);

				imgOut(x, ds1_ip->get_dim().sizes[1] - y - 1, z, 0) = col.r;
				imgOut(x, ds1_ip->get_dim().sizes[1] - y - 1, z, 1) = col.g;
				imgOut(x, ds1_ip->get_dim().sizes[1] - y - 1, z, 2) = col.b;
			}
		}
	}


	//Determine if one dataset is always larger
	int large_ds;
	if( 0 >= min && 0 <= max )
		large_ds = -1;
	else{

		layout loc;
		loc.arr_size = 3;
		loc.sizes = new int[3];
		loc.sizes[0] = 0;
		loc.sizes[1] = 0;
		loc.sizes[2] = 0;

		if(ds1_ip->get_indexed_point(loc).vals[var_ds1] 
			> ds2_ip->get_indexed_point(loc).vals[var_ds2]){

			large_ds = 1;
		}
		else
			large_ds = 2;
	}

	//Write border
	int border_start = ds1_ip->get_dim().sizes[0];

	for(int y=0; y < ds1_ip->get_dim().sizes[1]; y++){
		for(int x=border_start; x < border_start + border_width; x++){
			for(int z=0; z<ds1_ip->get_dim().sizes[2]; z++){

				imgOut( x, y, z, 0 ) = border_color.r;
				imgOut( x, y, z, 1 ) = border_color.g;
				imgOut( x, y, z, 2 ) = border_color.b;
			}
		}
	}

	//Write key
	int x_start = ds1_ip->get_dim().sizes[0] + border_width;

	long double k_height = ds1_ip->get_dim().sizes[1];
	long double unit = range / (k_height-1.0);

	for(int y=0; y<ds1_ip->get_dim().sizes[1]; y++){
		for(int x=x_start; x < x_start + key_width; x++){
			for(int z=0; z<ds1_ip->get_dim().sizes[2]; z++){
				rgb_color col = get_color_double_sided(min+unit*y, range, min, 0);

				imgOut( x, ds1_ip->get_dim().sizes[1] - y - 1, z, 0 ) = col.r;
				imgOut( x, ds1_ip->get_dim().sizes[1] - y - 1, z, 1 ) = col.g;
				imgOut( x, ds1_ip->get_dim().sizes[1] - y - 1, z, 2 ) = col.b;
			}
		}
	}

	//Construct file name
	string output_name(outprefix);
	output_name += "_difmap_wkey_";
	output_name += "var1_";
	output_name += itoa(var_ds1);
	output_name += "_var2_";
	output_name += itoa(var_ds2);

	//Write info to txt file
	output_info(output_name + ".txt", max, min, large_ds);

	//Write image
	imgOut.save_bmp((output_name + ".bmp").c_str());

	ret = output_name;
	return ret;
}



/* color map -------------------------------------------------------------*/

colmap::colmap(indexed_points *dataset1, indexed_points *dataset2, 
							int var_ds1, int var_ds2,
							ofstream &outfile, string outprefix) 
	: map(dataset1, dataset2, var_ds1, var_ds2, outfile, outprefix){
	
	return;
}

string colmap::get_res(){
	return compare();
}

string colmap::compare(){

	string ret("");

	indexed_points *ds1_ip = dynamic_cast<indexed_points*>(dataset1);
	indexed_points *ds2_ip = dynamic_cast<indexed_points*>(dataset2);

	if(ds1_ip == NULL || ds2_ip == NULL)
		return string("Improper Datasets");

	for(int j=1; j<=2; j++){
				
		long double max = get_max_value(j);
		long double min = get_min_value(j);

		long double range = max - min;

		int x_size = 0;
		int y_size = 0;
		int z_size = 0;

		if(j == 1){
			x_size = ds1_ip->get_dim().sizes[0];
			y_size = ds1_ip->get_dim().sizes[1];
			z_size = ds1_ip->get_dim().sizes[2];
		}	
		else{
			x_size = ds2_ip->get_dim().sizes[0];
			y_size = ds2_ip->get_dim().sizes[1];
			z_size = ds2_ip->get_dim().sizes[2];;
		}

		CImg<unsigned char> imgOut(x_size, y_size, z_size, 3, 0);

		for(int x=0; x<x_size; x++){
			for(int y=0; y<y_size; y++){
				for(int z=0; z<z_size; z++){

					layout loc;
					loc.arr_size = 3;
					loc.sizes = new int[3];
					loc.sizes[0] = x;
					loc.sizes[1] = y;
					loc.sizes[2] = z;

					long double pt_val;

					if(j == 1)
						pt_val = ds1_ip->get_indexed_point(loc).vals[var_ds1];
					else if(j == 2)
						pt_val = ds2_ip->get_indexed_point(loc).vals[var_ds2];

					rgb_color col = get_color_single_sided(pt_val, range, min);

					imgOut(x, y_size - y - 1, z, 0) = col.r;
					imgOut(x, y_size - y - 1, z, 1) = col.g;
					imgOut(x, y_size - y - 1, z, 2) = col.b;
				}
			}
		}

		//Construct file name
		string output_name(outprefix);
		output_name += "_colmap";

		if(j == 1){
			output_name += "_FIRST_";
			output_name += "var1_";
			output_name += itoa(var_ds1);	
		}
		else if(j == 2){
			output_name += "_SECOND_";
			output_name += "var2_";
			output_name += itoa(var_ds2);
		}

		output_name += ".bmp";

		ret += output_name + " ";

		imgOut.save_bmp(output_name.c_str());

	}

	return ret;

}

string colmap::get_metric_name(){
	return string("Color Map");
}

long double colmap::get_max_value(int ds){

	long double max = -std::numeric_limits<long double>::infinity();

	int size = dataset1->get_num_points();

	for(int i=0; i<size; i++){

		if(ds == 1){
			if(dataset1->get_point(i).vals[var_ds1] > max){
				max = dataset1->get_point(i).vals[var_ds1];
			}
		}

		if(ds == 2){
			if(dataset2->get_point(i).vals[var_ds2] > max){
				max = dataset2->get_point(i).vals[var_ds2];
			}
		}
	}

	return max;
}

long double colmap::get_min_value(int ds){

	long double min = std::numeric_limits<long double>::infinity();

	int size = dataset1->get_num_points();

	for(int i=0; i<size; i++){

		if(ds == 1){
			if(dataset1->get_point(i).vals[var_ds1] < min){
				min = dataset1->get_point(i).vals[var_ds1];
			}
		}

		if(ds == 2){
			if(dataset2->get_point(i).vals[var_ds2] < min){
				min = dataset2->get_point(i).vals[var_ds2];
			}
		}

	}

	return min;
}



/* sample correlation coefficient ----------------------------------------*/

scorco::scorco(points *dataset1, points *dataset2, int var_ds1, int var_ds2, 
				ofstream &outfile, string outprefix) 
				: metric(dataset1, dataset2, var_ds1, var_ds2, outfile, outprefix){
	
	return;
}

string scorco::get_res(){
	string ret;
	std::stringstream out;
	out.precision(10);
	out << compare();
	return out.str();
}

long double scorco::compare(){

	if(!check_card_match())
		return -std::numeric_limits<long double>::infinity();

	long double f1_ave = get_average_f1();
	long double f2_ave = get_average_f2();

	long double sum_f1 = 0;
	long double sum_sq_f1 = 0;

	long double sum_f2 = 0;
	long double sum_sq_f2 = 0;

	long double sum_mult = 0;

	for(int i=0; i<dataset1->get_num_points(); i++){
	
		long double f1 = (long double)dataset1->get_point(i).vals[var_ds1];
		long double f2 = (long double)dataset2->get_point(i).vals[var_ds2];
		
		if( !isnan(f1) && !isnan(f2) ){
			sum_f1 += f1;
			sum_sq_f1 += pow(f1, 2);
		
			sum_f2 += f2;
			sum_sq_f2 += pow(f2, 2);

			sum_mult += f1*f2;
		}

	}

	long double n = dataset1->get_num_points();

	long double res = 0;

	res = (n*sum_mult - (sum_f1*sum_f2)) / 
						( sqrt(n*sum_sq_f1 - pow(sum_f1,2)) * sqrt(n*sum_sq_f2 - pow(sum_f2,2)) );

	return res;

}

long double scorco::get_average_f1(){

	long double total = 0;	

	int num_indep = dataset1->get_num_indep_vars();

	for(int i=0; i<dataset1->get_num_points(); i++){
		total += (long double)dataset1->get_point(i).vals[num_indep];
	}

	return (total / ((long double)dataset1->get_num_points()));
}

long double scorco::get_average_f2(){
	
	long double total = 0;	

	int num_indep = dataset2->get_num_indep_vars();

	for(int i=0; i<dataset2->get_num_points(); i++){
		total += (long double)dataset2->get_point(i).vals[num_indep];
	}

	return (total / ((long double)dataset2->get_num_points()));
}


bool scorco::check_card_match(){	

	if(dataset1->get_num_points() == dataset2->get_num_points())
		return true;

	string err_ps("Mismatched number of coordinates: Dataset 1 = ");
	err_ps += itoa(dataset1->get_num_points());
	err_ps += ", Dataset 2 = ";
	err_ps += itoa(dataset2->get_num_points());

	error_out(err_ps);

	return false;
}

string scorco::get_metric_name(){
	return string("Sample Correlation Coefficient");
}



/* modelling efficiency --------------------------------------------------*/

modef::modef(points *dataset1, points *dataset2, int var_ds1, int var_ds2,
				ofstream &outfile, string outprefix) 
				: metric(dataset1, dataset2, var_ds1, var_ds2, outfile, outprefix){
	
	return;
}

string modef::get_res(){
	string ret;
	std::stringstream out;
	out.precision(10);
	out << compare();
	return out.str();
}

long double modef::compare(){

	if(!check_card_match())
		return -std::numeric_limits<long double>::infinity();

	long double f1_ave = get_average_f1();

	long double sum_sq_dif = 0;
	long double sum_sq_mean_dif = 0;

	for(int i=0; i<dataset1->get_num_points(); i++){

		long double f1 = (long double)dataset1->get_point(i).vals[var_ds1];
		long double f2 = (long double)dataset2->get_point(i).vals[var_ds2];
		
		if( !isnan(f1) && !isnan(f2) ){
			sum_sq_dif += pow(f1 - f2, 2);
			sum_sq_mean_dif += pow(f1 - f1_ave, 2);
		}

	}

	long double res = 0;

	res = 1 - ((sum_sq_dif) / (sum_sq_mean_dif));

	return res;

}

long double modef::get_average_f1(){

	long double total = 0;	

	int num_indep = dataset1->get_num_indep_vars();

	for(int i=0; i<dataset1->get_num_points(); i++){
		total += (long double)dataset1->get_point(i).vals[num_indep];
	}

	return (total / ((long double)dataset1->get_num_points()));
}

long double modef::get_average_f2(){
	
	long double total = 0;	

	int num_indep = dataset2->get_num_indep_vars();

	for(int i=0; i<dataset2->get_num_points(); i++){
		total += (long double)dataset2->get_point(i).vals[num_indep];
	}

	return (total / ((long double)dataset2->get_num_points()));
}


bool modef::check_card_match(){

	if(dataset1->get_num_points() == dataset2->get_num_points())
		return true;

	string err_ps("Mismatched number of coordinates: Dataset 1 = ");
	err_ps += itoa(dataset1->get_num_points());
	err_ps += ", Dataset 2 = ";
	err_ps += itoa(dataset2->get_num_points());

	error_out(err_ps);

	return false;
}

string modef::get_metric_name(){
	return string("Modelling Efficiency");
}




/* metric selection ------------------------------------------------------*/

void print_metrics(){
	
	int i=0;

	ofstream unset;
	string nullstr("");

	mse one(NULL, NULL, -1, -1, unset, nullstr);
	cout<<"\t"<<++i<<": "<<one.get_metric_name()<<"\n";

	rmse two(NULL, NULL, -1, -1, unset, nullstr);
	cout<<"\t"<<++i<<": "<<two.get_metric_name()<<"\n";

	scc three(NULL, NULL, -1, -1, unset, nullstr);
	cout<<"\t"<<++i<<": "<<three.get_metric_name()<<"\n";

	difmap four(NULL, NULL, -1, -1, unset, nullstr);
	cout<<"\t"<<++i<<": "<<four.get_metric_name()<<"\n";

	difmap_wkey five(NULL, NULL, -1, -1, unset, nullstr);
	cout<<"\t"<<++i<<": "<<five.get_metric_name()<<"\n";

	colmap six(NULL, NULL, -1, -1, unset, nullstr);
	cout<<"\t"<<++i<<": "<<six.get_metric_name()<<"\n";

	scorco seven(NULL, NULL, -1, -1, unset, nullstr);
	cout<<"\t"<<++i<<": "<<seven.get_metric_name()<<"\n";
	
	modef eight(NULL, NULL, -1, -1, unset, nullstr);
	cout<<"\t"<<++i<<": "<<eight.get_metric_name()<<"\n";

}

metric *metric_select(int m, points *ds1, points *ds2, int var_ds1, int var_ds2,
						ofstream &outfile, string outprefix){

	metric *ret;

	//Indexed points conversion
	// - For metrics that require indexed_points
	indexed_points *ds1_ip = dynamic_cast<indexed_points*>(ds1);
	indexed_points *ds2_ip = dynamic_cast<indexed_points*>(ds2);

	switch (m){

		case 1:
			ret = new mse(ds1, ds2, var_ds1, var_ds2, outfile, outprefix);
			break;

		case 2:
			ret = new rmse(ds1, ds2, var_ds1, var_ds2, outfile, outprefix);
			break;

		case 3:
			ret = new scc(ds1, ds2, var_ds1, var_ds2, outfile, outprefix);
			break;

		case 4:	
			ret = new difmap(ds1_ip, ds2_ip, var_ds1, var_ds2, outfile, outprefix);
			break;

		case 5:
			ret = new difmap_wkey(ds1_ip, ds2_ip, var_ds1, var_ds2, outfile, outprefix);
			break;

		case 6:
			ret = new colmap(ds1_ip, ds2_ip, var_ds1, var_ds2, outfile, outprefix);
			break;

		case 7:
			ret = new scorco(ds1, ds2, var_ds1, var_ds2, outfile, outprefix);
			break;

		case 8:
			ret = new modef(ds1, ds2, var_ds1, var_ds2, outfile, outprefix);
			break;

		default:
			ret = NULL;
			break;

	}

	return ret;

}

