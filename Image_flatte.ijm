function FlattenImage(imagename,type,molecule){
	run("Clear Results"); 	
	print(imagename);
	add_string = "Max Intensity";
	if (type == "SUM"){
		add_string = "Sum Slices";
	} else if(type == "MAX"){
		add_string = "Max Intensity"; 
	} else if(type == "MIN"){
		add_string = "Min Intensity";
	} else if(type == "AVG"){
		add_string = "Average Intensity"; 
	}else if(type == "STD"){
		add_string = "Standard Deviation"; 
	}else{
		add_string = "Median" ;
	}
	if(File.exists(folder+imagename)){
		open(folder+imagename);
		name=getTitle;
		run("Z Project...", "projection=["+add_string+"]");
		new_name = type+"_"+imagename;
		selectWindow(new_name);
		File.makeDirectory(folder+type);
		print("found file");
		save_path = folder+type+"/"+new_name;
		saveAs("tiff", save_path);
	}
	run("Close All");
}
//folder = "Gria1/";
root = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Anne-Sophie/";
molecule = "Gria2";
folder = root+molecule+"/";
list = getFileList(folder);
//type_list = newArray(["max","min","mean","median","sum","sd"]);
//type_list = [0,1,2,3,4,5]
num_type = 0;
if (num_type == 0){
		type = "SUM";
	} else if(num_type == 1){
		type = "MAX"; 
	} else if(num_type == 2){
		type = "MIN";
	} else if(num_type == 3){
		type = "AVG"; 
	}else if(num_type == 4){
		type = "STD" ;
	}else{
		type = "MED" ;
	}
//print(list[1]);
//Array.sort(list);
num_files = 0;
for (i = 0; i < list.length; i++){
	temp_fname =  split(list[i],".");
	if (temp_fname.length >= 2){
		if(temp_fname[1] == "lsm"){
			FlattenImage(list[i],type,molecule);
			print(list[i]);
			num_files++;
		}
	}
}
print(num_files);

