<?php
/// Copyright (C) by Wiktor `lakewik` Jezioro 1017
/// Chemical recation solver kernel module////

// Includes //

include "chemical_constants.php";
include "constants.php";
include "frontendcontent.php";
include "lang.php";

//// EXTERNAL LIBS //////////////////

include 'lib/finediff/finediff.php';

//////////////////////////////////////////////////////////////////////////////////////////////////////////

function find_atoms ($reaction) {
	// Variables:
	//		$recation - reaction or uncompleted reaction chain
	
	global $chemicalElements;
	
	/// Divide to atoms ///
	foreach ($chemicalElements as $value) {
		if (strpos($chemical_compound, $value) !== false) {
			$contains_iterator++;
			$substrats_atoms_array_parsed[$contains_iterator] = $value;
		}
	}
	
	return $substrats_atoms_array_parsed;
	
	
	
}

function char_at($str, $pos)
{
  return $str{$pos};
}

function divide_to_atoms ($chemical_compound, $use_acids_and_alcali, $options_array) {
	
	$tablica = Array();
	if ($use_acids_and_alcali) {
	
		if (contains_acid_moiety ($chemical_compound)) {
		
			$acid_moiety = get_acid_moiety($chemical_compound);
			$chemical_compound = str_replace($acid_moiety, "", $chemical_compound);
			//array_push ($tablica, $acid_moiety);
			$tablica[] = $acid_moiety;
		
		}
	
		if (strpos ($chemical_compound, "OH")) {// sprawdzamy czy zasada, jeśli tak wycinamy z analizy OH
	
			if ($options_array["apply_brackets"]) {
				
				if ( contains_brackets ($chemical_compound)) {
				
						$number_after_bracket = number_after_bracket ($chemical_compound);
						$str_replace_string = "(OH)".$number_after_bracket;
						
								
								
					}	else {
						$str_replace_string = "OH";
					}		
				
			} else {
				$str_replace_string = "OH";
			}
	
				$chemical_compound = str_replace($str_replace_string, "", $chemical_compound);
				array_push ($tablica, $str_replace_string);
	
		}
	}
	
	for ($j = 0; $j < strlen($chemical_compound); $j++)
		{
		$znak = char_at($chemical_compound, $j);
		if (ctype_upper($znak))
			{
			$atom = $znak;
			$j++;
			while ($j < strlen($chemical_compound))
				{
				$znak = char_at($chemical_compound, $j);
				
				if (ctype_lower($znak))
					{
						
					$atom .= $znak;
					$j++;
					}
				  else
					{
					$j--;
					break;
					}
				}
//echo $atom;
			array_push ($tablica, $atom);
			}
		}
		
		return $tablica;
	
	
}


function order_recation_products_array ($chemical_molecules_array) {
	
	$i == 0;
	foreach ($chemical_molecules_array as $value) {
					// if contains alcali

		if ( strpos($value, 'H2O') !== false OR strpos($value, 'H2')  !== false  OR strpos($value, '[not]H2[not]') !== false) { // sortujemy produkty bazując na wodzie
			//echo count ($chemical_molecules_array);
			$chemical_molecules_array_ready[0] = $value;
			if ($i == 0) {
				$firstiteration = true;
				//echo "test";
				//$reversed = array_reverse($input);
				$chemical_molecules_array_ready[count ($chemical_molecules_array)-1] = $chemical_molecules_array[1];
			}
			
			if ($i == 1) {
				//$chemical_molecules_array[0] = $chemical_molecules_array[1];
				$seconditeration = true;
			}
			
			break;
			
			
		}
		
	
		
		
		$i++;
	}
	
	if ($firstiteration) {
		
		 $chemical_molecules_array_ready = array_reverse( $chemical_molecules_array_ready);
		 return $chemical_molecules_array_ready;
	}
	
	if ($seconditeration) {
		return $chemical_molecules_array;
	}
	
	
	
	
}

function search_alone_gas_and_apply ($input_array) {
	// Funckja szukająca samotnuch gazów w tablicy i zastosowująca 2 atomy do samotnego gazu
	
	global $chemicalElementsProperty;
	
	$i = 0;
	foreach ($input_array as $value) {
		
		$value_ready[$i] = str_replace(" ", "", $value);
		
		if ($chemicalElementsProperty [$value_ready[$i]]["physicalState"] == "gas") {
			$value_ready[$i] = $value_ready[$i]."2"; // dodajemy żeby bylo 2 atomy
		}
		
		
		$i++;
	}
	
	
	return $value_ready;
	
	
}


function get_acid_from_acid_oxide ($acid_oxide) {
	
	global $acidOxideValue;
	$acid_oxide = remove_compounding_number ($acid_oxide);
	return $acidOxideValue[$acid_oxide];
	
}


function do_neutralization_reaction($substrats_array, $options_array, $substrats_types_array = null)
{ //   reakcja chemiczna między kwasem a zasadą lub tlenkiem metalu
	global $chemicalElements;
	global $FRONTEND;
	$i = 0;

	// / Do divide to atoms ///

	if (get_compound_type($substrats_array[1]) != "acid" AND get_compound_type($substrats_array[1]) != "acid_oxide") { ///  check and apply proper order for reaction
		$substrats_array_temp = $substrats_array;
		$substrats_array[1] = $substrats_array_temp[0];
		$substrats_array[0] = $substrats_array_temp[1];
	}

	foreach($substrats_array as $value) {
		if ($substrats_types_array == null) {
			$substrat_type[$i] = get_compound_type($value);
		}

		$substrat_atom_array[$i] = divide_to_atoms(change_acid_oxides_to_acids($value) , true, $options_array);

		// $substrats_parsed_array[$i] = change_acid_oxides_to_acids ($substrats_parsed_array[$i]);

		$substrats_parsed_array[$i] = order_compound_atoms_from_array(apply_valence($value, true));
		$substrats_parsed_array[$i] = change_acid_oxides_to_acids($substrats_parsed_array[$i]);

		//	echo $substrats_parsed_array[$i];

		$i++;
	}

	$reaction_chain = $substrats_array[0] . "+" . $substrats_array[1];

	// Do reaction
	// / 1

	if (get_neutralization_recation_solving_method($reaction_chain, $options_array) == 1) {

		// Do neutralization reaction with solution: ALCALI/HYDROXIDE + ACID -->> SALT + H2O

		$product_compounds_array[0] = $substrat_atom_array[0][0] . $substrat_atom_array[1][1];
		$product_compounds_array[1] = $substrat_atom_array[0][1] . $substrat_atom_array[1][0];
		$product_compounds_array = normalize_chemical_compounds_array($product_compounds_array);
		$product_compounds_array = order_recation_products_array($product_compounds_array);
	}

	// Do reaction
	// if (($substrat_type[0] == "metal" && $substrat_type[1] == "acid")
	//	OR (get_compound_type ($substrats_array[0]) == "acid" &&  get_compound_type ($substrats_array[1]) == "metal" )) {
	// / 2

	if (get_neutralization_recation_solving_method($reaction_chain, $options_array) == 2) {

		//	echo "test";
		// Do neutralization reaction with solution: METAL + ACID -->> SALT + H2

		$product_compounds_array[0] = $substrat_atom_array[0][0] . $substrat_atom_array[1][0];
		$product_compounds_array[1] = "H2"; // [not] - dont apply walence

		// $product_compounds_array = normalize_chemical_compounds_array ($product_compounds_array);
		// $product_compounds_array = order_recation_products_array($product_compounds_array);

	}

	// // 3

	if (get_neutralization_recation_solving_method($reaction_chain, $options_array) == 3) {

		//	echo "test";
		// Do neutralization reaction with solution: METAL OXIDE + ACID -->> SALT + H2O

		$product_compounds_array[0] = $substrat_atom_array[0][0] . $substrat_atom_array[1][0];
		$product_compounds_array[1] = $substrat_atom_array[0][1] . $substrat_atom_array[1][1];

		// $product_compounds_array = normalize_chemical_compounds_array ($product_compounds_array);
		// $product_compounds_array = order_recation_products_array($product_compounds_array);

	}

	// // 4

	if (get_neutralization_recation_solving_method($reaction_chain, $options_array) == 4) {
		$acid_oxides_method = true;

		// echo "test";
		// Do neutralization reaction with solution: ALCALI/HYDROXIDE + ACID OXIDE -->> SALT + H2O

		$acid = get_acid_from_acid_oxide($acid_oxide); // get acid fror acid oxide
		$product_compounds_array[0] = $substrat_atom_array[0][0] . $substrat_atom_array[1][1];
		$product_compounds_array[1] = $substrat_atom_array[0][1] . $substrat_atom_array[1][0];

		// $product_compounds_array = normalize_chemical_compounds_array ($product_compounds_array);
		// $product_compounds_array = order_recation_products_array($product_compounds_array);

		$acid_oxide_tmp = "";
	}

	// // 5

	if (get_neutralization_recation_solving_method($reaction_chain, $options_array) == 5) {
		$acid_oxides_method = true;
		$method5 = true;

		//	echo "test";
		// Do neutralization reaction with solution: METAL OXIDE + ACID OXIDE -->> SALT

		$product_compounds_array[0] = $substrat_atom_array[0][0] . $substrat_atom_array[1][0];

		//	$product_compounds_array[1] = $substrat_atom_array[0][1].$substrat_atom_array[1][1];
		// $product_compounds_array = normalize_chemical_compounds_array ($product_compounds_array);
		// $product_compounds_array = order_recation_products_array($product_compounds_array);

		$acid_oxide_tmp = "";
	}

	// / Parse products ///

	$i = 0;
	foreach($product_compounds_array as $value) {

		// print_r (apply_valence ( $value ));

		$product_compounds_array_parsed[$i] = order_compound_atoms_from_array(apply_valence($value, true)); /// DESIGN

		// $product_compounds_array_parsed[$i] = order_compound_atoms_from_array (  $value ); /// DESIGN
		// $product_compounds_array_parsed[$i] =  $value ; /// DESIGN
		// $product_compounds_array_parsed[$i] = apply_valence ( $value ); /// DESIGN

		$i++;
	}

	$product_compounds_array_parsed = normalize_chemical_compounds_array($product_compounds_array_parsed);
	$product_compounds_array_parsed = search_alone_gas_and_apply($product_compounds_array_parsed); //  Search alone Gas

	// return $substrat_atom_array;
	// return $substrat_atom_array;
	// apply_valence ($input)
	// return $substrats_parsed_array;
	// echo get_compound_type ($substrats_array[0]);
	// / Balancing reaction and return array of left and wight site

	if ($acid_oxides_method) {

		// / Check if method 3 or 4 ///
		// $acid_oxide_ready = search_acid_oxide_in_array ($substrats_array);
		// $acid_oxide_ready = search_acid_oxide_in_array ($substrats_array);

		if ($method5) {
			$substrats_metal_oxide_valence_tmp = apply_valence($substrats_array[0], true);
			$balanced_array = balance_reaction($substrats_metal_oxide_valence_tmp[0] . $substrats_metal_oxide_valence_tmp[1] . " + " . $substrats_array[1] . " = " . $product_compounds_array_parsed[0]);
		}
		else {
			$substrats_alcali_valence_tmp = apply_valence($substrats_array[0], true);
			$balanced_array = balance_reaction($substrats_alcali_valence_tmp[1] . $substrats_alcali_valence_tmp[0] . " + " . $substrats_array[1] . " = " . $product_compounds_array_parsed[0] . " + " . $product_compounds_array_parsed[1]);
		}
	}
	else {
		$balanced_array = balance_reaction($substrats_parsed_array[0] . " + " . $substrats_parsed_array[1] . " = " . $product_compounds_array_parsed[0] . " + " . $product_compounds_array_parsed[1]);
	}

	$balanced_array_without_html = preg_replace("/<.+>/sU", "", $balanced_array); /// Remove html tags

	// return $product_compounds_array_parsed;
	// / Exlode balanced to array
	// $balanced_array_exploded_left_site = explode ("+", $balanced_array[0]);
	// $balanced_array_exploded_right_site = explode ("+", $balanced_array[1]);
	// $balanced_array_exploded_with_arrows = insert_arrows ($balanced_array_exploded);  /// insert html arrows
	// $right_site = insert_arrows ($array)
	// print_r ($balanced_array_exploded_with_arrows);
	// print_r ($substrats_array);
	// //print_r ($product_compounds_array_parsed);
	// print_r ( $substrat_atom_array);
	// print_r ( $substrats_alcali_valence_tmp);
	// echo ($substrats_array[0]." + ".$substrats_array[1]." = ".$product_compounds_array_parsed[0]." + ".$product_compounds_array_parsed[1]);
	// / Molecule method enabled

	if ($options_array["molecule_enabled"]) {
		$output_array["molecule"] = $balanced_array;
	}

	// / Ions enabled

	if ($options_array["ion_enabled"]) {

		// $balanced_array_exploded = explode

		$options_array["in_water_valence"] = 1;
		$iterator = 0;
		$dissociated_array_0 = ion_dissociate($balanced_array[0], $type, $options_array);
		$dissociated_array_1 = ion_dissociate($balanced_array[1], $type, $options_array);
		foreach($dissociated_array_0 as $value) {
			$output_array["ion"].= $value;
			if ($iterator + 1 < count($dissociated_array_0)) {
				$output_array["ion"].= " + ";
			}

			$iterator++;
		}

		$output_array["ion"].= $FRONTEND["reaction_concat_sign"];
		$iterator = 0;
		foreach($dissociated_array_1 as $value) {
			$output_array["ion"].= $value;
			if ($iterator + 1 < count($dissociated_array_1)) {
				$output_array["ion"].= " + ";
			}

			$iterator++;
		}

		// $output_array["ion"][] = ion_dissociate ($balanced_array[1], $type, $options_array);

	}

	// / Shorten Ions enabled

	if ($options_array["ion_cut_enabled"]) {

		// $balanced_array_exploded = explode

		$options_array["in_water_valence"] = 1;
		$iterator = 0;
		$dissociated_array_0 = ion_dissociate($balanced_array[0], $type, $options_array);
		$dissociated_array_1 = ion_dissociate($balanced_array[1], $type, $options_array);
		foreach($dissociated_array_0 as $value) {
			$output_array["ion_cut"].= $value;
			if ($iterator + 1 < count($dissociated_array_0)) {
				$output_array["ion_cut"].= " + ";
			}

			$iterator++;
		}

		$output_array["ion_cut"].= $FRONTEND["reaction_concat_sign"];
		$iterator = 0;
		foreach($dissociated_array_1 as $value) {
			$output_array["ion_cut"].= $value;
			if ($iterator + 1 < count($dissociated_array_1)) {
				$output_array["ion_cut"].= " + ";
			}

			$iterator++;
		}

		$output_array["ion_cut"] = shorten_ions($output_array["ion_cut"], $options_array_ions_shorten);
	}

	$output_array["errors"] = null;
	return $output_array;
}


function do_slurry_reaction ($input_array, $options_array) {
	
	/*
	Method 1 - SALT1 + SALT2 = SALT3 + SALT4
	Method 2 - SOULABLE ACID + SOULABLE SALT = OTHER SOULABLE ACID + OTHER SALT PRECIPITATE
	Method 3a - ALCALI (GROUP 2) + SALT = OTHER ALCALI + OTHER SALT PRECIPITATE
	Method 3b - ALCALI (GROUP 1 and 2) + SALT = OTHER HYDROXIDE PRECIPITATE + OTHER SALT 
	
	*/
	
	$substrats_array = $input_array["substrats"];
	
	$i = 0;
	foreach ($substrats_array as $value) {
		
			
		if ($substrats_types_array == null) {
			$substrat_type[$i] = get_compound_type ($value);
		}
		
		$substrat_atom_array[$i] = divide_to_atoms ($value, true, $options_array);
		
	
		//$substrats_parsed_array[$i] = change_acid_oxides_to_acids ($substrats_parsed_array[$i]);
		
		
		
		$substrats_parsed_array[$i] =  apply_valence($value, true); 
		//$substrats_parsed_array[$i] = change_acid_oxides_to_acids ($substrats_parsed_array[$i]);
		
	//	echo $substrats_parsed_array[$i];
		
		$i++;
	}
	
	$reaction_chain = $substrats_array[0]."+".$substrats_array[1];
	
	
	
	if ($input_array["reaction_type"] == "auto") {
		/// Try auto-obtain slurry reaction type ///
		
		
	} elseif ($input_array["reaction_type"] == "1") {
		$reaction_type = 1;
	} elseif ($input_array["reaction_type"] == "2") {
		$reaction_type = 2;
	} elseif ($input_array["reaction_type"] == "3a") {
		$reaction_type = 3;
	} elseif ($input_array["reaction_type"] == "3b") {
		$reaction_type = 4;
	} elseif ($input_array["reaction_type"] == "4") { //// SPECIAL - This reaction can be solved using two methods - 3 and 4
		$reaction_type = 5;
	}
	
	
	if ($reaction_type == 1) {
		// 	Method 1 - SALT1 + SALT2 = SALT3 + SALT4
		
			$product_compounds_array[0] = $substrat_atom_array[1][1].$substrat_atom_array[0][0];
			$product_compounds_array[1] = $substrat_atom_array[0][1].$substrat_atom_array[1][0];
			$product_compounds_array = normalize_chemical_compounds_array ($product_compounds_array);
		
		
	}
	
	if ($reaction_type == 2) {
		// 	Method 2 - SOULABLE ACID + SOULABLE SALT = OTHER SOULABLE ACID + OTHER SALT PRECIPITATE
		
			$product_compounds_array[0] = $substrat_atom_array[1][1].$substrat_atom_array[0][0];
			$product_compounds_array[1] = $substrat_atom_array[0][1].$substrat_atom_array[1][0];
			$product_compounds_array = normalize_chemical_compounds_array ($product_compounds_array);
		
		
	}
	
	
	if ($reaction_type == 3) {
		// 	Method 3a - ALCALI (GROUP 2) + SALT = OTHER ALCALI + OTHER SALT PRECIPITATE
		
			$product_compounds_array[0] = $substrat_atom_array[1][1].$substrat_atom_array[0][0];
			$product_compounds_array[1] = $substrat_atom_array[0][1].$substrat_atom_array[1][0];
			$product_compounds_array = normalize_chemical_compounds_array ($product_compounds_array);
		
		
	}
	
	if ($reaction_type == 4) {
		// 	Method 3b - ALCALI (GROUP 1 and 2) + SALT = OTHER HYDROXIDE PRECIPITATE + OTHER SALT 
		
			$product_compounds_array[0] = $substrat_atom_array[1][1].$substrat_atom_array[0][0];
			$product_compounds_array[1] = $substrat_atom_array[0][1].$substrat_atom_array[1][0];
			$product_compounds_array = normalize_chemical_compounds_array ($product_compounds_array);
		
		
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if ($reaction_type == 5) {  ////// DEPRECATED ////
		// 	Method 3a and 3b
		// 	Method 3a - ALCALI (GROUP 2) + SALT = OTHER ALCALI + OTHER SALT PRECIPITATE
		// 	Method 3b - ALCALI (GROUP 1 and 2) + SALT = OTHER HYDROXIDE PRECIPITATE + OTHER SALT 
		
		
			$product_compounds_array[0] = $substrat_atom_array[1][1].$substrat_atom_array[0][0];
			$product_compounds_array[1] = $substrat_atom_array[0][1].$substrat_atom_array[1][0];
			$product_compounds_array = normalize_chemical_compounds_array ($product_compounds_array);
			
			
			$product_compounds_array[0] = $substrat_atom_array[1][1].$substrat_atom_array[0][0];
			$product_compounds_array[1] = $substrat_atom_array[0][1].$substrat_atom_array[1][0];
			$product_compounds_array = normalize_chemical_compounds_array ($product_compounds_array);
		
		
	}
	
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	$i = 0;
	foreach ($product_compounds_array as $value) {
		$product_compounds_array_parsed[$i] = order_compound_atoms_from_array (  apply_valence ( $value, true)); /// DESIGN
		$i++;
	}
	
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
	$product_compounds_array_parsed = normalize_chemical_compounds_array ($product_compounds_array_parsed);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  	$balanced_array = balance_reaction($substrats_parsed_array[0][1].$substrats_parsed_array[0][0] . " + " . $substrats_parsed_array[1][1].$substrats_parsed_array[1][0] . " = " . $product_compounds_array_parsed[0] . " + " . $product_compounds_array_parsed[1]);
  
	//print_r ($product_compounds_array_parsed);
	print_r ($substrats_parsed_array);
	print_r ( $balanced_array);
	
}

function get_error_number_description ($error, $options_array) {
	
	if ($error == "E31") {
		$error_message = "W zapisie reakcji został podany metak, który jest mało aktywny, w związku z czym taka reakcja nie może zachodzić.";
		if ($options_array["display_solution"] == 1) {
			$solution_message = "Należy użyć matelu aktywnego. Metale aktywne to: ";
		}
	}
	
}

function why_cant_do_reaction ($substrats_array, $options_array) {
	
	
}


function get_slurry_recation_solving_method ($substrats_array, $options_array) {
	
	if (get_compound_type ($substrats_array[0]) == "salt" AND get_compound_type ($substrats_array[1]) == "salt") {
		if (if_compound_soluble_in_water ($substrats_array[0]) AND if_compound_soluble_in_water ($substrats_array[1])) {
			return 1;
		}
		
	} elseif (get_compound_type ($substrats_array[0]) == "salt" AND get_compound_type ($substrats_array[1]) == "acid") {
		if (if_compound_soluble_in_water ($substrats_array[0])) {
			return 2;
		}
		
	} elseif (get_compound_type ($substrats_array[0]) == "acid" AND get_compound_type ($substrats_array[1]) == "salt") {
		if (if_compound_soluble_in_water ($substrats_array[1])) {
			return 2;
		}
		
	} elseif (get_compound_type ($substrats_array[0]) == "alcali" AND get_compound_type ($substrats_array[1]) == "salt") {
		
		$atoms_array = divide_to_atoms ($substrats_array[0], true, $options_array);
		if (get_chemical_element_property ("group", $atoms_array[1]) == 2) {
			return 3;
		} elseif (get_chemical_element_property ("group", $atoms_array[1]) == 2 || get_chemical_element_property ("group", $atoms_array[1]) == 1) {
			return 4;
		}
		
	} elseif (get_compound_type ($substrats_array[0]) == "salt" AND get_compound_type ($substrats_array[1]) == "alcali") {
		
		$atoms_array = divide_to_atoms ($substrats_array[1], true, $options_array);
		if (get_chemical_element_property ("group", $atoms_array[1]) == 2) {
			return 3;
		} elseif (get_chemical_element_property ("group", $atoms_array[1]) == 2 || get_chemical_element_property ("group", $atoms_array[1]) == 1) {
			return 4;
		}
		
		
	}
	
	//print_r ($substrats_array);
	
	
	
}

function get_slurry_reaction_compounds_types_chain ($substrats_array, $options_array) {
		global $FRONTEND;
	

	
	/// 1 
	if (get_slurry_recation_solving_method ($substrats_array, $options_array) == 1) {	
			return "Sól 1 Rozpuszczalna + Sól 2 Rozpuszczalna = Sól 3 Trudno rozpuszczalna (Osad) + Sól 4 Rozpuszczalna (Przeszącz)";
	}
	
	/// 2 
	if (get_slurry_recation_solving_method ($substrats_array, $options_array) == 2) {	
			return "Kwas rozpuszczalny + Sól rozpuszczalna = Inny kwas rozpuszczalny + Inna sól (osad)";
	}
	
	/// 3 
	if (get_slurry_recation_solving_method ($substrats_array, $options_array) == 3) {	
			return "Zasada (Grupa 1) + Sól = Inna zasada + Inna sól (osad) &darr;";
	}
	
	/// 4
	if (get_slurry_recation_solving_method ($substrats_array, $options_array) == 4) {	
			return "Zasada (Grupa 1 i 2) + Sól = Inny wodorotlenek (osad) &darr; + Inna sól";
	}


	
}


function is_reaction_ready_to_make () {
	
	
	return true;
}

function get_neutralization_reaction_compounds_types_chain ($reaction_chain, $options_array) {
		global $FRONTEND;
	

	
	/// 1 
 if (get_neutralization_recation_solving_method ($reaction_chain, $options_array) == 1) {
			
		// Do neutralization reaction with solution: ALCALI/HYDROXIDE + ACID -->> SALT + H2O

        return "Zasada/Wodorotlenek + Kwas = Sól + Woda";
		
		}

	//Do reaction 
/// 2	
 if (get_neutralization_recation_solving_method ($reaction_chain, $options_array) == 2) {
        return "Metal aktywny + Kwas = Sól + Wodór";

		
		}	
//// 3		
 if (get_neutralization_recation_solving_method ($reaction_chain, $options_array) == 3) {
        return "Tlenek metalu + Kwas = Sól + Woda";

		
		}	
		
//// 4		
 if (get_neutralization_recation_solving_method ($reaction_chain, $options_array) == 4) {
			        return "Zasada/Wodorotlenek + Tlenek kwasowy/Bezwodnik kwasowy = Sól + Woda";


		}	
//// 5		
 if (get_neutralization_recation_solving_method ($reaction_chain, $options_array) == 5) {
			        return "Tlenek metalu + lenek kwasowy/Bezwodnik kwasowy = Sól";

		
		}	
		
		

		
		//return $substrats_type_array[0].$substrats_type_array[1].$products_type_array[0].$products_type_array[1];
		//return $substrats_names_ready.$FRONTEND["reaction_concat_sign"].$products_names_ready;
	
}


function get_neutralization_recation_solving_method ($reaction_chain, $options_array) {
	
	
	$reaction_chain = str_replace (" ", "", $reaction_chain);
		$reaction_chain_exploded = explode("=", $reaction_chain);
		$substrats_parsed_array = explode("+", $reaction_chain_exploded[0]);
	//$solving_method = null;
	
	/// 1
	// Do neutralization reaction with solution: ALCALI/HYDROXIDE + ACID -->> SALT + H2O
	 if (is_hydroxide($substrats_parsed_array[0]) && get_compound_type ($substrats_parsed_array[1]) == "acid")
		  {
			return 1;
			
		} elseif (get_compound_type ($substrats_parsed_array[0]) == "acid" && is_hydroxide($substrats_parsed_array[1])) {
						return 1;
		} 
		
		
		/// 2
	// Do neutralization reaction with solution: ACTIVE METAL + ACID -->> SALT + H2
		
		elseif (get_compound_type ($substrats_parsed_array[0]) == "acid" && get_compound_type ($substrats_parsed_array[1]) == "metal") {
						return 2;
		} elseif (get_compound_type ($substrats_parsed_array[1]) == "acid" && get_compound_type ($substrats_parsed_array[0]) == "metal") {
						return 2;
		}
		
			/// 3
	// Do neutralization reaction with solution: METAL OXIDE + ACID -->> SALT + H2O
		
		elseif (get_compound_type ($substrats_parsed_array[0]) == "acid" && get_compound_type ($substrats_parsed_array[1]) == "metal_oxide") {
						return 3;
		} elseif (get_compound_type ($substrats_parsed_array[1]) == "acid" && get_compound_type ($substrats_parsed_array[0]) == "metal_oxide") {
						return 3;
		}
		
		
		/// 4
	// Do neutralization reaction with solution: ALCALI/HYDROXIDE + ACID OXIDE -->> SALT + H2O
	
	
	 elseif (is_hydroxide($substrats_parsed_array[0]) && get_compound_type ($substrats_parsed_array[1]) == "acid_oxide")
		  {
			return 4;
			
		} elseif (get_compound_type ($substrats_parsed_array[0]) == "acid_oxide" && is_hydroxide($substrats_parsed_array[1])) {
						return 4;
		} 
	
	
	
	/// 5
	// Do neutralization reaction with solution: METAL OXIDE + ACID OXIDE -->> SALT
	
		
		elseif (get_compound_type ($substrats_parsed_array[0]) == "acid_oxide" && get_compound_type ($substrats_parsed_array[1]) == "metal_oxide") {
						return 5;
		} elseif (get_compound_type ($substrats_parsed_array[1]) == "acid_oxide" && get_compound_type ($substrats_parsed_array[0]) == "metal_oxide") {
						return 5;
		} else {
			return false;
		}
		
	

	
	
	
	
	
	//return $solving_method;
	
	
}


function get_any_reaction_compounds_types_chain ($reaction_chain, $options_array) {
		global $FRONTEND;
	
		$i = 0;
		
		$substrats_names_ready = null;
		$products_names_ready = null;
		
		$reaction_chain = str_replace (" ", "", $reaction_chain);
		$reaction_chain_exploded = explode("=", $reaction_chain);
		$substrats_exploded_array = explode("+", $reaction_chain_exploded[0]);
		$products_exploded_array = explode("+", $reaction_chain_exploded[1]);
	
	foreach ($substrats_exploded_array as $value) {
		$substrats_type_array[$i] = change_compounds_name_to_user_friendly_names (get_compound_type ($value), $options_array);
		$substrats_names_ready .= $substrats_type_array[$i]." ";
	
		$i++;
		//echo $substrats_type_array[$i-1];
		if (count($substrats_exploded_array) > $i && $value !== "") {
			$substrats_names_ready .= "&plus; ";
		}
	}
	
	$i = 0;
	
	foreach ($products_exploded_array as $value) {
		$products_type_array[$i] = change_compounds_name_to_user_friendly_names (get_compound_type ($value), $options_array);
		$products_names_ready .= $products_type_array[$i]." ";
		$i++;
		
		if (count($products_exploded_array) > $i && $value !== "") {
			$products_names_ready .= "&plus; ";
		}
		
	}
	

		
		//return $substrats_type_array[0].$substrats_type_array[1].$products_type_array[0].$products_type_array[1];
		return $substrats_names_ready.$FRONTEND["reaction_concat_sign"].$products_names_ready;
	
}


function change_acid_oxides_to_acids ($molecules, $options_array = null) {
	
		
	if (is_array($molecules)) {
		
	} else {
		
		if (contains_acid_oxide ($molecules)) {
		
		$ready = get_acid_from_acid_oxide ($molecules);
			return $ready;
		} else {
			return $molecules;
		}
	}
	
	
}


function contains_acid_oxide ($molecule, $options_array = null) {
	
	//global $acidOxideValueReverse;
	global $acidOxideValue;
	$molecule = remove_compounding_number ($molecule);
	
	//echo $molecule;
	
	//if (array_search($molecule,$acidOxideValueReverse) !== false) {
	if (array_key_exists($molecule, $acidOxideValue)) {
		return true;
	} else {
		return false;
	}
		
	
}

function search_acid_oxide_in_array ($array, $options_array = null) {
	
	
}


function remove_compounding_number ($chemical_compound) {
	// Usuwanie liczby cząsteczkowej z poczętku cząsteczki
	$i = 0;
	$chemical_compound_parsed = strip_tags($chemical_compound, "<sub><sup>"); /// Remove html tags
	$chemical_compound_parsed = preg_replace('/\s+/', ' ', trim($chemical_compound_parsed));
	$chemical_compound_parsed = str_replace(" ", "", $chemical_compound_parsed); /// Remove whitespaces
	
	while (is_numeric($chemical_compound_parsed[0+$i])) {
		$chemical_compound_parsed[0+$i] = "";
		$i++;
	}
	
	return $chemical_compound_parsed;
}

function get_compounding_number ($chemical_compound) {
	// Zawraca liczbę cząsteczkową
	$i = 0;
	$chemical_compound_parsed = strip_tags($chemical_compound);
	while (is_numeric($chemical_compound_parsed[0+$i])) {
		$compounding_number_tmp .= $chemical_compound_parsed[0+$i];
		$i++;
	}
	
	return $compounding_number_tmp;
}


function insert_arrows ($array) {
	// Dodawanie strzałek ulatniania sie gazów lub osiadania osadu
	
	$i = 0;
	foreach ($array as $value) {
		
		$value_parsed_tmp = remove_compounding_number ($value);
		//echo $value_parsed_tmp;
		if ((string)$value_parsed_tmp == 'H<sub>2</sub>') {
			echo "dobrze";
			$compounding_number = get_compounding_number ($value); // Restore deleted compounding number
			$final_array[$i] = $compounding_number.$value_parsed_tmp."&#8593;"; // insert HTML ASCI Upper Arrow
		} else {
			echo "zle";
			$final_array[$i] = $value;
		}
		
		
		$i++;
	}
	
	return $final_array;
	
}

function get_reaction_type ($reaction_chain, $options_array) {
	
	
}


function strSlice($str, $start, $end) {
    $end = $end - $start;
    return substr($str, $start, $end);
}

function balance_reaction ($reaction_chain, $options_array = null) {
	
	$pathToRhinoJar = '/var/www/html/reactionsolver/engine/rhino1.7.7.1/lib/rhino-1.7.7.1.jar';
	//$javascriptFile = '/var/www/html/reactionsolver/engine/chemicalbalance.js';
	$javascriptFile = '/var/www/html/reactionsolver/engine/chemicalbalance_phantom.js';
	//$output = shell_exec("java -jar $pathToRhinoJar $javascriptFile ".'"'.$reaction_chain.'"');
	$output = shell_exec("QT_QPA_PLATFORM=offscreen phantomjs $javascriptFile ".'"'.$reaction_chain.'"');
   
   //print_r ($output);
   
    $output_array_1 = explode ("#$#", $output);
    $output_array = explode ("<=>", $output_array_1[1]);
   
	return $output_array;

	
	
}


function parse_double_hydrogens ($chemical_compound) {
	
	$chemical_compound_parsed = str_replace ("OHH", "H2O", $chemical_compound);
	return $chemical_compound_parsed;
	
}

function normalize_h2o ($chemical_compound) {
	
	if($chemical_compound != strip_tags($chemical_compound)) {
		// contains HTML
		$chemical_compound = strip_tags ($chemical_compound);
		$chemical_compound_parsed = str_replace ("OH2", "H<sub>2</sub>O", $chemical_compound);
	} else {
		$chemical_compound_parsed = str_replace ("OH2", "H2O", $chemical_compound);
	}
	
	//echo $chemical_compound ;
	//echo "test2";
	return $chemical_compound_parsed;
	
}

function normalize_chemical_compounds_array ($array) {
	
	foreach ($array as &$value) {
		$value = parse_double_hydrogens ($value);
		$value = normalize_h2o ($value);
		
	}
	//print_r ($array);
	return $array;
	
}



function contains_hydrogen_at_begin ($chemical_compound) {
	if ($chemical_compound[0] == "H") {
		
		return true;
		
	} else {
		
		return false;
		
	}
	
}

function parse_acid_moiety_ion_to_html ($scid_moiety_ion) {
	
	
}

function is_hydroxide ($chemical_compound) {
	
	if (strpos($chemical_compound, 'OH') !== false) {
		return true;
	} else {
		return false;
	}
	
}

function get_compound_type ($chemical_compound) {
	
	//echo $chemical_compound;
	
	if (contains_acid_moiety ($chemical_compound)) {
		if (contains_hydrogen_at_begin ($chemical_compound)) {
					$is_acid = true;
					return "acid";
					
		}
		
	} 
	
	if (contains_acid_oxide ($chemical_compound) && $is_acid != true) {
			return "acid_oxide";
	}
	
	
	
	
	
	if (strpos($chemical_compound, 'OH') !== false) {
		if (if_compound_soluble_in_water ($chemical_compound)) {
			return "alcali";
		} else {
		
			return "hydroxide";
		
		}
		
	} 
	
	
	
	
	
	if (strpos($chemical_compound, "H2O") !== false) {
		
		return "water";
		
	}
	
	
	if (contains_acid_moiety ($chemical_compound)) {
		if (contains_metal_at_begin ($chemical_compound)) {
					return "salt";
		}
		
	}
	
	
	if (get_chemical_element_property ("type", $chemical_compound) == "metal") {
		return "metal";
	}
	
	
	
	if (strpos($chemical_compound, 'O') !== false) { /// check if contains OXYGEN - this is optimization
	
		$atoms_array = divide_to_atoms ($chemical_compound, false, null);
		
		if (get_chemical_element_property ("type", $atoms_array[0]) == "metal" && $atoms_array[1] == "O") {
			return "metal_oxide";
		}
		
		
	}
	
	
	
	
}





function change_compounds_name_to_user_friendly_names ($chemical_compound, $options_array = null) {
	
	global $LANG;
	
	$chemical_compound = strtolower($chemical_compound);
	
	return $LANG[$chemical_compound];
	
}


function contains_metal_at_begin ($chemical_compound) {
	global $chemicalElements;
	
	/*
	
	foreach ($chemicalElements as $value) {
	
		if ($chemical_compound[0] == $value) {
			if (get_chemical_element_property ("type", $chemical_compound[0]) == "metal") {
				return true;
				break;
			}
		} elseif ($chemical_compound[0].$chemical_compound[1] == $value) {
			if (get_chemical_element_property ("type", $chemical_compound[0].$chemical_compound[1]) == "metal") {
				return true;
				break;
			}
			
		}
		
	}

	*/	
	
	if (get_chemical_element_property ("type", $chemical_compound[0]) == "metal" AND ctype_lower($chemical_compound[1]) !== true) {
				return true;
		} elseif (get_chemical_element_property ("type", $chemical_compound[0].$chemical_compound[1]) == "metal") {
				return true;
		}
		
	
}


function get_minimal_valence ($atom) {
	
	
}

function do_dissociation ($chemical_compound) {
	
	
	$dissociated = str_replace ("OH", "OH-", $chemical_compound);
	
	return $dissociated;
	
	
	
}

function is_single_atom ($input) {
	
	/// Check if input is single atom
	
}


function if_compound_soluble_in_water ($chemical_compound) {
	
	global $solubilityArray;
	
	if (strpos($chemical_compound, 'OH') !== false) { // check if hydroxide
		
		
		
		$truncated_chemical_compound = trim(preg_replace('/\s*\([^)]*\)/', '', $chemical_compound)); // Remove dual or more hydroxide groups in brackets
		 
		$truncated_chemical_compound = str_replace ("OH", "", $truncated_chemical_compound); // Remole dingle hydroxides groups
		
		$truncated_chemical_compound = str_replace (" ", "", $truncated_chemical_compound); // Remove whitespaces
		
		$truncated_chemical_compound = preg_replace('/[0-9]+/', '', $truncated_chemical_compound); // Remove all numbers
		
		if ($solubilityArray ["OH"] [$truncated_chemical_compound] == "1") {
			
			return true;
			
		} else {
			
			return false;
			
		}
		
		
	}
	
	///////////// FOR ANY OTHER ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	$atoms_array = divide_to_atoms (remove_compounding_number($chemical_compound), true, $options_array);
	
	if ($solubilityArray [$atoms_array[0]] [$atoms_array[1]] == "1") {
			
			return true;
			
		} else {
			
			return false;
			
		}
	
	
	
	//echo $truncated_chemical_compound;
	
	
	
}

function add_brackets ($input) {
	
	return "(".$input.")";
	
}


function contains_oxide ($chemical_compound) {
	
	$chemical_compound_exploded = explode ("+", $chemical_compound);
	
	$i = 0;
	foreach ($chemical_compound_exploded as $value) {
		
		$value_tmp = remove_compounding_number($value);
		
		if (strpos(get_compound_type($value_tmp), "oxide") !== false) {
			
			return true;
			break;
			
		}
		
		$i++;
	}
	
}


function get_oxides_in_chemical_compound ( $chemical_compound, $options_array = null) {
	
	$chemical_compound_exploded = explode ("+", $chemical_compound);
	
	$i = 0;
	foreach ($chemical_compound_exploded as $value) {
		
		$value_tmp = remove_compounding_number($value);
		//echo $value_tmp;
		if (strpos(get_compound_type($value_tmp), "acid_oxide") !== false OR strpos(get_compound_type($value_tmp), "metal_oxide") !== false) {
			
			if ($options_array["get_oxides"] == "one") {
				$oxide = $value;
				return $oxide;
				break;
			} elseif ($options_array["get_oxides"] == "all") {
				$oxides_array[] = $value;
			}
			
		}
		
		$i++;
	}
	
	if ($options_array["get_oxides"] == "all") {
	
		return $oxides_array;
		
	}
	
	
}

function shorten_ions ($ions_equation_chain, $options_array) {
	global $FRONTEND;
	
	
	///// Function for cutting the same ions from left and right site in ions equation chain ////////
	
	//$ions_equation_chain_clear = str_replace(" ", "", $ions_equation_chain);
	
	//$ions_equation_chain = implode('+',array_unique(explode('+', $ions_equation_chain)));
		
	$ions_equation_chain_exploded = explode("&rarr;", $ions_equation_chain);
	
	/*
	
*/
//echo $ions_equation_chain;

//$from_text = "3OH<sup>&#43;</sup> + PO4<sup>3&#43;</sup> + 3H<sup>&#45;</sup>";
//$to_text = "PO4<sup>3&#43;</sup> + H2O";

$from_text = $ions_equation_chain_exploded[0];
$to_text = $ions_equation_chain_exploded[1];

$opcodes = FineDiff::getDiffOpcodes($from_text, $to_text, FineDiff::wordDelimiters);

FineDiff::renderCommonsFromOpcodes($from_text, $opcodes);

$shorte_array = FineDiff::$commons;

//print_r(FineDiff::$commons);

$ions_equation_chain_exploded[0] = str_replace ($shorte_array, "", $ions_equation_chain_exploded[0]);
$ions_equation_chain_exploded[1] = str_replace ($shorte_array, "", $ions_equation_chain_exploded[1]);
$ions_equation_chain_exploded = array_map('trim', $ions_equation_chain_exploded); /// trin whitespaces
//print_r( $ions_equation_chain_exploded);
$ions_equation_chain_exploded = str_replace ("   ", " ", $ions_equation_chain_exploded);
$ions_equation_chain_exploded = str_replace ("  ", " ", $ions_equation_chain_exploded);
$ions_equation_chain_exploded = str_replace (" ", " + ", $ions_equation_chain_exploded);
$ions_equation_chain_exploded = str_replace (" + + + ", " + ", $ions_equation_chain_exploded);
$ions_equation_chain_exploded = str_replace (" + + ", " + ", $ions_equation_chain_exploded);

return $ions_equation_chain_exploded[0].$FRONTEND["reaction_concat_sign"].$ions_equation_chain_exploded[1];



	
}

function order_atoms_for_dissociation ($atoms_array, $options_array = null) {
	
	
		if (contains_acid_moiety($atoms_array[0])) { ///  check and apply proper order for reaction 
		
		$atoms_array_temp = $atoms_array;
		
		$atoms_array[1] = $atoms_array_temp[0];
		$atoms_array[0] = $atoms_array_temp[1];
	 
	
	} elseif (is_hydroxide($atoms_array[0])) {
		
		$atoms_array_temp = $atoms_array;
		
		$atoms_array[1] = $atoms_array_temp[0];
		$atoms_array[0] = $atoms_array_temp[1];
		
		
	}
	
	return $atoms_array;
	
}

function get_atom_number_in_compound ($chemical_compound, $options_array) {
			
			$chemical_compound = str_replace(" ", "", $chemical_compound);
	$chemical_compound = strip_tags($chemical_compound); // Delete HTML tags
	$this_options_array["get_oxides"] = "one";
						$value3 = $chemical_compound;
						if ($options_array["divide_atoms_include_special_compounds"]) { /// jeżeli reszty kwasowe lub końcowki OH sa włączone do rozdziału
							$divide_to_atoms_tmp1 = false;
						} else {
							$divide_to_atoms_tmp1 = true;
						}
			$atoms_array = divide_to_atoms($value3, $divide_to_atoms_tmp1, $options_array);
			foreach($atoms_array as $value) {
			
			
			
			if (strpos($value3, $value) !== false) { // określamy pozycję atomu w cząsteczce
				$pos = strpos($value3, $value);
				if (strlen($value) == 2) { //// chemical element name length check
					$multiplier = 2;
				}
				elseif (strlen($value) == 1) {
					$multiplier = 1;
				}

				$atoms_number_tmp[$i] = "";

				// / checking if any special molecule

				if (contains_brackets($value3) && $value == "OH") {
					$atoms_number_tmp[$i] = number_after_bracket($value3);

					// echo "nawiasy";

				}
				elseif (contains_acid_moiety($value) && contains_brackets($value3)) { // check if probably hydroxide or alcali
					$atoms_number_tmp[$i] = number_after_bracket($value3);

					//	echo "reszta kwasowa";

				}
				elseif (contains_acid_moiety($value)) { // check if probably hydroxide or alcali
					$atoms_number_tmp[$i] = 1;

					// echo "reszta kwasowa";

				}
				else {
					$i2 = 0;
					while (is_numeric($value3[$pos + $multiplier + $i2])) {
						$atoms_number_tmp[$i].= $value3[$pos + $multiplier + $i2];
						$i2++;
					}
				}
			}
			if ($atoms_number_tmp[$i] == "") {
				$atoms_number_tmp[$i] = 1;
			}
			$final_atoms_number_array[$value] = $atoms_number_tmp[$i];

			$i++;
			
			//echo $i;
		}
		
	return  ($final_atoms_number_array);
}

 
function ion_dissociate($chemical_compound, $type, $options_array = null)
{
	
	$chemical_compound = str_replace(" ", "", $chemical_compound);
	$chemical_compound = strip_tags($chemical_compound); // Delete HTML tags
	$this_options_array["get_oxides"] = "one";
	
	///// Check if contains acid oxide
	if (contains_oxide($chemical_compound)) {
		$not_dissociating_array[0] = get_oxides_in_chemical_compound($chemical_compound, $this_options_array);
		$chemical_compound = str_replace($not_dissociating_array[0], "", $chemical_compound);
	} 
	
	

	$i = 0;
	$i3 = 0;

	// echo $not_dissociating_array[0];
	// / !!!!! /// // explode chemical compound

	$chemical_compound_array = explode("+", $chemical_compound);
	$chemical_compound_array = array_filter($chemical_compound_array);
	//print_r ($chemical_compound_array);
	foreach($chemical_compound_array as $value3) {
//echo "kotek";
		// //////////
//echo $value3;
	//	print_r($chemical_compound_array);
	
		$atoms_array = order_atoms_for_dissociation (divide_to_atoms($value3, true, $options_array));
	
		if (strpos($value3, "H") !== false && count ($atoms_array) == 1 ) {
		 
		 /// contains only hydrogen - add to not-dissociating array ////
		 
		 $not_dissociating_array[] = $value3;
		 
		 $chemical_compound = str_replace($value3, "", $chemical_compound);
		 
		 $value3 = "";
		 $atoms_array = order_atoms_for_dissociation (divide_to_atoms($value3, true, $options_array));
		} 
		
		//// Check if contains water
		elseif (contains_water ( remove_compounding_number($value3) )) {
		//echo "kotek";
		$not_dissociating_array[] = $value3;

		
		
		//$not_dissociating_array[1] = "H2O";

				 
		//$chemical_compound = str_replace($not_dissociating_array[1], "", $chemical_compound);
		$chemical_compound = str_replace($value3, "", $chemical_compound);
		$value3 = "";
		$atoms_array = null;
		}
		
		else {
	
	
		$atoms_array = order_atoms_for_dissociation (divide_to_atoms($value3, true, $options_array));
		$compounding_number = get_compounding_number($value3);
		//echo $compounding_number ;
		if ($compounding_number == "") {
			$compounding_number = 1;
		}

		// $options_array = Array();


		// / numbers to iona ex. XNO3

		foreach($atoms_array as $value) {
			
			
			
			if (strpos($value3, $value) !== false) { // określamy pozycję atomu w cząsteczce
				$pos = strpos($value3, $value);
				if (strlen($value) == 2) { //// chemical element name length check
					$multiplier = 2;
				}
				elseif (strlen($value) == 1) {
					$multiplier = 1;
				}

				$atoms_number_tmp[$i] = "";

				// / checking if any special molecule

				if (contains_brackets($value3) && $value == "OH") {
					$atoms_number_tmp[$i] = number_after_bracket($value3);

					// echo "nawiasy";

				}
				elseif (contains_acid_moiety($value) && contains_brackets($value3)) { // check if probably hydroxide or alcali
					$atoms_number_tmp[$i] = number_after_bracket($value3);

					//	echo "reszta kwasowa";

				}
				elseif (contains_acid_moiety($value)) { // check if probably hydroxide or alcali
					$atoms_number_tmp[$i] = 1;

					// echo "reszta kwasowa";

				}
				else {
					$i2 = 0;
					while (is_numeric($value3[$pos + $multiplier + $i2])) {
						$atoms_number_tmp[$i].= $value3[$pos + $multiplier + $i2];
						$i2++;
					}
				}
			}

			$i++;
			//echo $i;
		}
		
		}

		// $i = 0;

		foreach($atoms_array as $value) {
			if ($i3 % 2 == 0) { /// determine charge sign
				$znak = "&#43;"; /// plus
			}
			else {
				$znak = "&#45;"; /// minus
			}

			if ($atoms_number_tmp[$i3] == 0 OR $atoms_number_tmp[$i3] == "") {
				$atoms_number_tmp[$i3] = 1;
			}

			$final_compounding_number[$i3] = $compounding_number * $atoms_number_tmp[$i3]; // set number at compound begin ex. XHNO3
			if ($final_compounding_number[$i3] == 1) {
				$final_compounding_number[$i3] = "";
			}

			if (contains_brackets($value3) && $value == "OH") {
				$molecule_valence = 1;
				if ($molecule_valence == 1) {
					$molecule_valence = "";
				}

				$ions_array[$i3] = $final_compounding_number[$i3] . "" . $value . "<sup>" . $molecule_valence . $znak . "</sup>";
			}
			elseif (contains_acid_moiety($value) && contains_brackets($value3)) { // check if probably hydroxide or alcali
				$molecule_valence = get_acid_moiety_valence($value);
				if ($molecule_valence == 1) {
					$molecule_valence = "";
				}

				$ions_array[$i3] = $final_compounding_number[$i3] . "" . $value . "<sup>" . $molecule_valence . $znak . "</sup>";
			}
			elseif (contains_acid_moiety($value)) { // check if probably hydroxide or alcali
				$molecule_valence = get_acid_moiety_valence($value);
				if ($molecule_valence == 1) {
					$molecule_valence = "";
				}

				$ions_array[$i3] = $final_compounding_number[$i3] . $value . "<sup>" . $molecule_valence . $znak . "</sup>";
			}
			else {
				if ($options_array["in_water_valence"] == "1") {
					$molecule_valence = get_chemical_element_property("valence_water", $value);
				}
				else {
					$molecule_valence = get_chemical_element_property("valence", $value);
				}

				if ($molecule_valence == 1) {
					$molecule_valence = "";
				}

				$ions_array[$i3] = $final_compounding_number[$i3] . $value . "<sup>" . $molecule_valence . $znak . "</sup>";
			}

			$i3++;
		}
	}
	
	//print_r ($atoms_number_tmp);
	//print_r ($final_compounding_number);
//	print_r ($atoms_array);

	///// Append not dissociating array ////////
	foreach($not_dissociating_array as $value) {
		$ions_array[] = $value;
	}
	

	// print_r ($ions_array);

	return (array_filter($ions_array));
	if ($type = "full") {
	}
	elseif ($type = "reduced") {
	}
}

function count_atoms_in_molecule () {
	
	
	
}

function contains_brackets ($input) {
	
	if (strpos($input, '(') !== false AND strpos($input, ')') !== false) { 
		return true;
	} else {
		return false;
	}
	
}

function order_compound_atoms_from_array ($chemical_atoms_array) {  // ustawianie w odpowiedniej kolejnosci atomow w związku chemicznym
	
	foreach ($chemical_atoms_array as $value) {
		if (strpos($value, 'OH') !== false) {  // if contains alcali
			
			//if (strpos($value, '(') !== false OR strpos($value, ')') !== false) { 
			
				//$nawias_begin = strpos($value, "(");
				//$nawias_end = strpos($value, ")");
	
	
				//$ile_razy_wziete = $value[$nawias_end+1];
				//$atend =  "".$value."".$ile_razy_wziete;
			
			//} else {
				$atend = $value;
			//}
			
			
			if ($i == 0) {
				$atbegin = $chemical_atoms_array[1];
			}
			
			if ($i == 1) {
				$atbegin = $chemical_atoms_array[0];
			}
			
			break;
			
			
		} elseif (contains_acid_moiety ($value)) {
			
			if (strpos($value, '(') !== false OR strpos($value, ')') !== false) { 
			
				$nawias_begin = strpos($value, "(");
				$nawias_end = strpos($value, ")");
	
	
				$ile_razy_wziete = $value[$nawias_end+1];
				$atend =  "(".get_acid_moiety ($value).")".$ile_razy_wziete;
			
			} else {
				$atend =  get_acid_moiety ($value);
			}
			
			if ($i == 0) {
				$atbegin = $chemical_atoms_array[1];
			}
			
			if ($i == 1) {
				$atbegin = $chemical_atoms_array[0];
			}
			
			break;
			
		} elseif (array_contains_water($value)) {
			echo "woda";
			$atbegin = "H2";
			$atend = "O";
		} else {
			
			$atbegin = $chemical_atoms_array[0].$chemical_atoms_array[1];
			
		}
		
		
		$i++;
	}
	
	//print_r ($chemical_atoms_array);
	
	return $atbegin.$atend;
	
	
}

function array_contains_water ($compound_array) {
	
	if ( $compound_array[0] == "H2" AND $compound_array[1] == "O"  ) {
		return true;
	} elseif ( $compound_array[1] == "H2" AND $compound_array[0] == "O" ) {
		return true;
	} else  {
		return false;
	}
	
}


function contains_water ($chemical_compound) {
	
	//echo $chemical_compound;
	
	if (strpos($chemical_compound, "H2O") !== FALSE) {
		return true;
	} else  {
		return false;
	}
	
}

function is_non_oxyacid ($chemical_compound) {
	
	
}

function is_non_oxysalt ($chemical_compound) {
	
	if (contains_acid_moiety($chemical_compound) != true) {
		return true;
	} else {
		return false;
	}
	
}


function apply_valence ($input, $with_water = false, $options_array = null) {
	/// Funkcja nanosząca i zastosowująca wartościowości
	
	//$input2 = $input;
	$i = 0;
	$atoms_array = divide_to_atoms ($input, true, $options_array);
	/// Iterate divided atoms array
	foreach ($atoms_array as $value) {
		/// Check if it is acid moiety or other non-atom
		if (contains_acid_moiety ($value)) {
			$atoms_valence_array[$i] = get_acid_moiety_valence($value);
		} elseif ($value == "OH") {// sprawdzamy czy wodorotlenek, jeśli tak wycinamy z analizy OH
	
			$atoms_valence_array[$i] = 1;
	
		} else {
			if ($with_water) {
				$atoms_valence_array[$i] = get_chemical_element_property("valence_water", $value);
			} else {	
				$atoms_valence_array[$i] = get_chemical_element_property("valence", $value);
			}
		}

		$i++;
	}
	
	// Count least common multiple 
	$least_common_multiple = lcd($atoms_valence_array,1); 
	
	/// Iterate atoms valence array
	$i = 0;
	foreach ($atoms_valence_array as $value) {
		
		$atoms_count_tmp = $least_common_multiple/$value;
		if ($atoms_count_tmp < 1) {
			$atoms_count_tmp = "";
		}
		
		
		/// Check if it is acid moiety or other non-atom
		if (contains_acid_moiety ($atoms_array[$i])) {
			if ($least_common_multiple/$value > 1) {
				
				//echo "test";
				
				$molecule_ready[$i] = add_brackets ($atoms_array[$i]);
				$molecule_ready[$i] .= $least_common_multiple/$value;
				
			} 
			
			else {
				
				$molecule_ready[$i] = $atoms_array[$i];
				//$molecule_ready[$i] .= $least_common_multiple/$value;
				
			}
		} elseif ($atoms_array[$i] == "OH") {// sprawdzamy czy zasada, jeśli tak wycinamy z analizy OH
	
			if ($least_common_multiple/1 > 1) {
				
				$molecule_ready[$i] = add_brackets ($atoms_array[$i]);
				$molecule_ready[$i] .= $least_common_multiple/$value;
				
			} 
			
			else {
				
				$molecule_ready[$i] = $atoms_array[$i];
				//$molecule_ready[$i] .= $least_common_multiple/1;
				
			}
			
	
		}
		
		else {
			$atoms_count_tmp = $least_common_multiple/$value;
			if ($atoms_count_tmp > 1) {
				$molecule_ready[$i] = $atoms_array[$i].$atoms_count_tmp;
			} else { // nic nie stoi to znaczy 1
				$molecule_ready[$i] = $atoms_array[$i];	
			}
		}

		$i++;
		//echo $atoms_count_tmp;
	}
	
	//echo $least_common_multiple;
	//print_r ($atoms_valence_array);
		//echo $input;


	//if (strpos($input2, "[not]")) {
	//	echo "test";
	//	return $input2;
	//} else {
		return $molecule_ready;
	//}
	//return $atoms_valence_array;
	//return $atoms_array;
	
	
}


function get_acid_moiety_valence($acid_moiety) {
	global $acidMoietyValence;
	return $acidMoietyValence[$acid_moiety];
	
}


function get_chemical_element_property ($property, $chemical_element) {
	
	global $chemicalElementsProperty;
	
	
	return $chemicalElementsProperty[$chemical_element][$property];

	
	
}


function contains_acid_moiety ($chemical_compound) {
	// Sprawdza, czy zawiera resztę kwasową
	global $acidMoiety;
	
	if (array_contains_part_of_string ($chemical_compound, $acidMoiety)) {
		return true;
	} else {
		return false;
		
	}
	
	
}


function get_acid_moiety ($chemical_compound) {
	// Zwraca resztę kwasową
	global $acidMoiety;
	$acidMoietyReady =  get_array_containing_part_of_string ($chemical_compound, $acidMoiety);
	return $acidMoietyReady;
	
}


function array_contains_part_of_string ($string, $array) {	
// check if string contains some words from array
	foreach ($array as $value) {
		//if (strstr($string, $url)) { // mine version
		if (strpos($string, $value) !== FALSE) { // Yoshi version
			return true;
			break;
			
		}
	}
  return false;
	
}



function get_array_containing_part_of_string ($string, $array) {	
// check if string contains some words from array
	foreach ($array as $value) {
		//if (strstr($string, $url)) { // mine version
		if (strpos($string, $value) !== FALSE) { // Yoshi version
			return $value;
			break;
			
		}
	}
  return false;
	
}



function contains($str, array $arr)
{
    foreach($arr as $a) {
        if (stripos($str,$a) !== false) return true;
    }
    return false;
}

//  Greatest Common Factor 
function gcf($a, $b) { 
	return ( $b == 0 ) ? ($a):( gcf($b, $a % $b) ); 
}


/// Least Common Multiple of two numbers
function lcm($a, $b) { 
	return ( $a / gcf($a,$b) ) * $b; 
}


///  Least Common Multiple of more than two numbers
function lcd($ar) {
               
        if (count($ar) > 1) {
			$ar[] = lcm( array_shift($ar) , array_shift($ar) );
			return lcd( $ar );
		} else {
			return $ar[0];
		}
}



function get_compound_types_chain ($reaction_chain) {
	
	$reaction_chain = str_replace (" ", "", $reaction_chain); 
	$exploded_chain =  explode ("=", $substrats_string_parsed);
	
	$substrats = $exploded_chain[0];
	$products = $exploded_chain[1];
	
	
	$substrats_exploded = explode ("+", $substrats);
	$products_exploded = explode ("+", $products);
	
	$i = 0;
	foreach ($substrats_exploded as $value) {
		
		$substrat_type[$i] = get_compound_type ($value);
		$i++;
	}
	
	
	$i = 0;
	foreach ($products_exploded as $value) {
		
		$product_type[$i] = get_compound_type ($value);
		$i++;
	}
	
	$substrat_type_polish = parse_to_polish_compounds_type_names ($substrat_type);
	$product_type_polish = parse_to_polish_compounds_type_names ($product_type);
	
}



function parse_to_polish_compounds_type_names ($compounds) {
	
	if (is_array($compounds)) {
		
	} else {
		
		
	}
	
}

function number_after_bracket ($chemical_compound) {
	
	$nawias_begin = strpos($chemical_compound, "(");
	$nawias_end = strpos($chemical_compound, ")");
	
	$ile_razy_wziete = $chemical_compound[$nawias_end+1];
	
	return $ile_razy_wziete;
	
}

function number_after_bracket_advanced ($chemical_compound, $molecule) {
	
	
}


function count_atoms_in_chemical_compound ($chemical_compound) {
	
	
	$nawias_begin = strpos($chemical_compound, "(");
	$nawias_end = strpos($chemical_compound, ")");
	
	$roznica = $nawias_end-$nawias_begin-1;
	
	$ile_razy_wziete = $chemical_compound[$nawias_end+1]-1;
	
	
	
	
	$i = 0;
	$numberPos = 0;
	while ($i < $ile_razy_wziete) {
		$i2 = 0;
		$i++;
		while ($i2 < $roznica) {
			$i2++;
			$chemical_compound[$nawias_end+$numberPos] = $chemical_compound[$nawias_end-$i2];
			
			$numberPos++;
		}
		
		
	}
	
	//$nawias_end+1
	//strlen($chemical_compound)
	
	//echo $chemical_compound;
	//echo $nawias_begin." ".$nawias_end;
	
	$atomsNumber = countAtoms($chemical_compound);
	if ($chemical_compound[0] == "" OR is_numeric($chemical_compound[0]) == false) {
		$compoundsNumber = 1;
	} else {
		$compoundsNumber = $chemical_compound[0];
	}
	

	
	$atomsNumberArray = $atomsNumber[0]*$compoundsNumber;
	
	return $atomsNumberArray;
	
}
$text = "2Fe(OH)3";


function countAtoms($text)
	{
	$tablica = Array();
	for ($j = 0; $j < count($tablicaAtomow); $j++)
		{
		$tablica[$j] = 0;
		}

	for ($j = 0; $j < strlen($text); $j++)
		{
		$znak = char_at($text, $j);
		if (ctype_upper($znak))
			{
			$atom = $znak;
			$number = "";
			$j++;
			while ($j < strlen($text))
				{
				$znak = char_at($text, $j);
				if (ctype_lower($znak))
					{
					$atom+= $znak;
					$j++;
					}
				  else
					{
					break;
					}
				}

			if (is_numeric(char_at($text, $j)))
				{
				while ($j < strlen($text))
					{
					$znak = char_at($text, $j);
					if (is_numeric($znak))
						{
						$number+= $znak;
						$j++;
						}
					  else
						{
						$j--;
						break;
						}
					}
				}
			  else
				{
				$j--;
				$number = "1";
				}

			//$index = array_search($atom, $tablicaAtomow);
			$tablica[0]+= $number;
			}
		}

	return $tablica;
	}


	
	
	

//print_r (   count_atoms_in_chemical_compound ($text));



///// CP TASKS - ZADANIA NA STĘŻENIE PROCENTOWE ///////

function do_cp_task ($input_array, $options_array) {
	
	
	
	
}

//// PRAWO ZACHOWANIA MASY ////

function atom_mass_task ($input_array, $options_array) {
	
	global $CONST_atomicMass;
	
	if ($input_array["task_type"] == 1) {
		/// Oblicz % zawartość pierwiastków w związkach chemicznych
		
		//// Liczymy poszczególne atomy
		
		
	}
	
	
	if ($input_array["task_type"] == 2) {
		/// Znając masę atomu [pierwiastek] która wynosi [liczba w notacji wykładniczej] oblicz jego masę atomową
		
		//// Liczymy proporcję
		/*
		
		1u - 0,166 * 10 -23
		x - 0,167 * 10 -23
		
		0,167 * 10 -23 = 0,166 * 10 -23 x
		
		
		
		*/
		
		// Calculate formula ///
		
		$atomic_mass = $input_array["atom_mass"]/$CONST_atomicMass;
		
		return $atomic_mass;
		
	}
	
	if ($input_array["task_type"] == 2) {
		/// Znając masę atomową [pierwiastek] która wynosi [liczba w unitach] oblicz jego masę atomu
		
		//// Liczymy proporcję
		/*
		
		*/
		
		// Calculate formula ///
		
		$atom_mass = $input_array["atomic_mass"]*$CONST_atomicMass;
		
		return $atom_mass;
		
	}
	
	
	
}

////////////////////////////////////////////////////////////////////////////////////////

function do_task_from_template ($input_array) {
	
	if ($input_array["template_id"] == "1") {
		
		$ms = $input_array["ms"];
		$mr = $ms + $input_array["mh2o"];
		
		$score = round ( ($ms * 100)/$mr,  $input_array["round_dec_num"]);
		
		$frac1 = '<div class="frac">
    <span>ms x 100%</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">mr</span>
    
</div>';
		
		$calculations = '	<hr><h2>Dane:</h2>
		Masa substancji rozpuszczonej (ms): <b>'.$ms.'g</b>
						<br>Masa roztworu (mr): <b>'.$ms.'g + '.$input_array["mh2o"].'g = '.$mr.'</b>
						<br>Stężenie roztworu (Cp): <b>?</b>
												<hr><h2>Szukane:</h2>
						<br>Stężenie roztworu (Cp): <b>?</b>

						<hr><h2>Wzór:</h2>
						<br> <b>Cp = &nbsp;'.$frac1.'</b>
						<hr><h2>Obliczenia:</h2>
						<br>Cp = <div class="frac">
    <span>'.$ms.'g &bull;	 100%</span>
    <span class="symbol">/</span>
    <span class="bottom">'.$mr.'g</span>
</div> = <b>'.$score.'%</b>';
		
	}
	
	
	if ($input_array["template_id"] == "2") {
		
		$cp = $input_array["cp"];
		$mr = $input_array["mr"];
		
		$score = round ( ($cp * $mr)/100,  $input_array["round_dec_num"]);
		
		$frac1 = '<div class="frac">
    <span>ms &bull; 100%</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">mr</span>
    
</div>';
		
		$calculations = '	<hr><h2>Dane:</h2>
						<br>Masa roztworu (mr): <b>'.$mr.'</b>
						<br>Stężenie roztworu (Cp): <b>'.$cp.'%</b>
						<hr><h2>Szukane:</h2>
								Masa substancji rozpuszczonej (ms): <b>?</b>

						<hr><h2>Wzory:</h2>
						<br> <b>Cp = &nbsp;'.$frac1.'</b>
						<br>Wzór po przekształceniu:&nbsp;&nbsp;
						<b>ms = 
						<div class="frac">
    <span>Cp &bull; mr</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">100%</span>
    
</div></b>

						<hr><h2>Obliczenia:</h2>
						<br>ms = <div class="frac">
    <span>'.$cp.'% &bull;	'.$mr.'g</span>
    <span class="symbol">/</span>
    <span class="bottom">100%</span>
</div> = <b>'.$score.'g</b>';
		
	}
	
	
	if ($input_array["template_id"] == "3") {
		
		$cp1 = $input_array["cp"];
		$mr = $input_array["mr"];
		$ms_added = $input_array["ms_added"];
		
		$ms = round ( ($cp1 * $mr)/100,  $input_array["round_dec_num"]);
		
		$ms_final = $ms_added + $ms;
		$mr_final = $ms_added + $mr;
		
		$score = round ( ($ms_final * 100)/$mr_final,  $input_array["round_dec_num"]);
		
		
		
		$frac1 = '<div class="frac">
    <span>ms &bull; 100%</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">mr</span>
    
</div>';
		
		$calculations = '	<hr><h2>Dane:</h2>
						<br>Masa roztworu (mr): <b>'.$mr.'</b>
						<br>Ilość dodanej substancji: <b>'.$ms_added.'</b>
						<br>Pierwotne stężenie roztworu (Cp<sub>1</sub>): <b>'.$cp1.'%</b>
						<hr><h2>Szukane:</h2>
								<br>Stężenie roztworu po dodaniu substancji (Cp<sub>2</sub>): <b>?</b>

						<hr><h2>Wzory:</h2>
						<br> <b>Cp = &nbsp;'.$frac1.'</b>
						<br>Wzór po przekształceniu:&nbsp;&nbsp;
						<b>ms = 
						<div class="frac">
    <span>Cp &bull; mr</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">100%</span>
    
</div></b>

						<hr><h2>Obliczenia:</h2>
						<br><u>Obliczamy pierwotną masę substancji (ms):</u>
						<br>ms = <div class="frac">
    <span>'.$cp1.'% &bull;	'.$mr.'g</span>
    <span class="symbol">/</span>
    <span class="bottom">100%</span>
</div> = <b>'.$ms.'g</b>
<br><br>
<u>Dodajemy do obliczonej masy substancji masę substancji dodanej do roztworu:</u>
<br>'.$ms.'g + '.$ms_added.'g = <b>'.$ms_final.'g</b><br>
<br>
<u>Dodajemy do masy roztworu masę substancji dodanej do roztworu:</u>
<br>'.$mr.'g + '.$ms_added.'g = <b>'.$mr_final.'g</b><br>
<br>
<u>Liczymy finalne stężenie procentowe:</u>
<br>Cp<sub>2</sub> = <div class="frac">
    <span>'.$ms_final.'g &bull;	 100%</span>
    <span class="symbol">/</span>
    <span class="bottom">'.$mr_final.'g</span>
</div> = <b>'.$score.'%</b>

';
		
	}
	
	if ($input_array["template_id"] == "4") {
		
		$cp1 = $input_array["cp"];
		$mr = $input_array["mr"];
		$ms_added = $input_array["ms_added"];
		
		$ms = round ( ($cp1 * $mr)/100,  $input_array["round_dec_num"]);
		
		$ms_final = $ms_added + $ms;
		$mr_final = $ms_added + $mr;
		
		$score = round ( ($ms_final * 100)/$mr_final,  $input_array["round_dec_num"]);
		
		
		
		$frac1 = '<div class="frac">
    <span>ms &bull; 100%</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">mr</span>
    
</div>';
		
		$calculations = '	<hr><h2>Dane:</h2>
						<br>Masa roztworu (mr): <b>'.$mr.'</b>
						<br>Ilość dodanej substancji: <b>'.$ms_added.'</b>
						<br>Pierwotne stężenie roztworu (Cp<sub>1</sub>): <b>'.$cp1.'%</b>
						<hr><h2>Szukane:</h2>
								<br>Stężenie roztworu po dodaniu substancji (Cp<sub>2</sub>): <b>?</b>

						<hr><h2>Wzory:</h2>
						<br> <b>Cp = &nbsp;'.$frac1.'</b>
						<br>Wzór po przekształceniu:&nbsp;&nbsp;
						<b>ms = 
						<div class="frac">
    <span>Cp &bull; mr</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">100%</span>
    
</div></b>

						<hr><h2>Obliczenia:</h2>
						<br><u>Obliczamy pierwotną masę substancji (ms):</u>
						<br>ms = <div class="frac">
    <span>'.$cp1.'% &bull;	'.$mr.'g</span>
    <span class="symbol">/</span>
    <span class="bottom">100%</span>
</div> = <b>'.$ms.'g</b>
<br><br>
<u>Dodajemy do obliczonej masy substancji masę substancji dodanej do roztworu:</u>
<br>'.$ms.'g + '.$ms_added.'g = <b>'.$ms_final.'g</b><br>
<br>
<u>Dodajemy do masy roztworu masę substancji dodanej do roztworu:</u>
<br>'.$mr.'g + '.$ms_added.'g = <b>'.$mr_final.'g</b><br>
<br>
<u>Liczymy finalne stężenie procentowe:</u>
<br>Cp<sub>2</sub> = <div class="frac">
    <span>'.$ms_final.'g &bull;	 100%</span>
    <span class="symbol">/</span>
    <span class="bottom">'.$mr_final.'g</span>
</div> = <b>'.$score.'%</b>

';
		
	}
	
	
	if ($input_array["template_id"] == "5") {
		
		$cp1 = $input_array["cp1"];
		$cp2 = $input_array["cp2"];
		$mr1 = $input_array["mr1"];
		$mr2 = $input_array["mr2"];
		$water_added = $input_array["water_added"];
		
		$ms1 = round ( ($cp1 * $mr1)/100,  $input_array["round_dec_num"]);
		$ms2 = round ( ($cp2 * $mr2)/100,  $input_array["round_dec_num"]);
		
		$ms3 = $ms1 + $ms2;
		
		$mr3 = $mr1 + $mr2 + $water_added;
		
		$score = round ( ($ms3 * 100)/$mr3,  $input_array["round_dec_num"]);
		
		
		
		
		
		$calculations = '	<hr><h2>Dane:</h2>
						<br>Gęstość wodyy : <b>1 g/cm<sup>3</sup></b> (w obliczeniach przyjmujemy gęstość wody która wynosi 1 g/cm<sup>3</sup> w temperaturze 4 C)
						<br>Masa pierwszego roztworu (mr<sub>1</sub>): <b>'.$mr1.'g</b>
						<br>Masa drugiego roztworu (mr<sub>2</sub>): <b>'.$mr2.'g</b>
						<br>Ilość dolanej wody: <b>'.$water_added.'cm<sup>3</sup></b>
						<br>Pierwotne stężenie procentowe pierwszego roztworu (Cp<sub>1</sub>): <b>'.$cp1.'%</b>
						<br>Pierwotne stężenie procentowe drugiego roztworu (Cp<sub>1</sub>): <b>'.$cp2.'%</b>
				
						<hr><h2>Szukane:</h2>
								<br>Nowe stężenie procentowe roztworu po zmieszaniu dwóch roztworów (Cp<sub>3</sub>): <b>?</b>

						<hr><h2>Wzory:</h2>
						<br> <b>Cp = &nbsp;'.$frac1.'</b>
						<br>Wzór po przekształceniu:&nbsp;&nbsp;
						<b>ms = 
						<div class="frac">
    <span>Cp &bull; mr</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">100%</span>
    
</div></b>

<br>Wzór na gęstość:
	<b>d = 
						<div class="frac">
    <span>m</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">V</span>
    
</div></b>						<br>Wzór po przekształceniu: (będzie nam potrzebny do obliczenia masy wody):&nbsp;&nbsp;
	<b>m = d &bull; V</b>


						<hr><h2>Obliczenia:</h2>
						<br>1. <u>Obliczamy pierwotną masę substancji w pierwszym roztworze (ms):</u>
						<br>ms = <div class="frac">
    <span>'.$cp1.'% &bull;	'.$mr1.'g</span>
    <span class="symbol">/</span>
    <span class="bottom">100%</span>
</div> = <b>'.$ms1.'g</b>
	<br>2. <u>Obliczamy pierwotną masę substancji w drugim roztworze (ms):</u>
						<br>ms = <div class="frac">
    <span>'.$cp2.'% &bull;	'.$mr2.'g</span>
    <span class="symbol">/</span>
    <span class="bottom">100%</span>
</div> = <b>'.$ms2.'g</b>
<br><br>
3. <u>Dodajemy obliczone masy dwóch substancji:</u>
<br>'.$ms1.'g + '.$ms2.'g = <b>'.$ms3.'g</b><br>
<br>
4. <u>Obliczamy masę wody podstawiając dane z wcześniej przekształconego wzoru na gęstość:</u>
<br>1 g/cm<sup>3</sup> &bull; '.$water_added.'cm<sup>3</sup> = <b>'.$water_added.'g</b><br>

<br>
5. <u>Dodajemy do masy roztworu pierwszego masę roztworu drugiego i masę wody:</u>
<br>'.$mr1.'g + '.$mr2.'g  + '.$water_added.'g  = <b>'.$mr3.'g</b><br>
<br>
6. <u>Liczymy finalne stężenie procentowe otrzymanego roztworu:</u>
<br>Cp<sub>3</sub> = <div class="frac">
    <span>'.$ms3.'g &bull; 100%</span>
    <span class="symbol">/</span>
    <span class="bottom">'.$mr3.'g</span>
</div> = <b>'.$score.'%</b>

';
		
	}
	
	
	if ($input_array["template_id"] == "6") {
		
		$ms = $input_array["ms"];
		$mr = $ms + $input_array["mh2o"];
		
		$score = round ( ($ms * 100)/$mr,  $input_array["round_dec_num"]);
		
		$frac1 = '<div class="frac">
    <span>ms x 100%</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">mr</span>
    
</div>';
		
		$calculations = '	<hr><h2>Dane:</h2>
		Masa substancji rozpuszczonej (ms): <b>'.$ms.'g</b>
						<br>Masa roztworu (mr): <b>'.$ms.'g + '.$input_array["mh2o"].'g = '.$mr.'</b>
						<br>Stężenie roztworu (Cp): <b>?</b>
												<hr><h2>Szukane:</h2>
						<br>Stężenie roztworu (Cp): <b>?</b>

						<hr><h2>Wzór:</h2>
						<br> <b>Cp = &nbsp;'.$frac1.'</b>
						<hr><h2>Obliczenia:</h2>
						<br>Cp = <div class="frac">
    <span>'.$ms.'g &bull;	 100%</span>
    <span class="symbol">/</span>
    <span class="bottom">'.$mr.'g</span>
</div> = <b>'.$score.'%</b>';
		
	}
	
	
	if ($input_array["template_id"] == "12") {
		
		///// USED VARIABLES ////////////////////////////////
		
		$weight = $input_array["input_weight_1"];
		$element_1 = $input_array["input_element_1"];
		$element_2 = $input_array["input_element_2"];
		
		///////////////////////////////////////////////
		
		$element_1_parsed = $element_1;
		$element_2_parsed = $element_2;
		
		$atomic_mass_1 =  round (calculate_molar_mass ($element_1_parsed, $options_array = null),  $input_array["round_dec_num_atomic_mass"]);
		$atomic_mass_2 =  round (calculate_molar_mass ($element_2_parsed, $options_array = null),  $input_array["round_dec_num_atomic_mass"]);
		
		$atoms_element_1 = round ( ($weight * MOL)/$atomic_mass_1,  $input_array["round_dec_num"]);
		$score = round ( ($atomic_mass_2 * $atoms_element_1)/MOL,  $input_array["round_dec_num"]);
		
		$frac1 = '<div class="frac">
    <span>ms x 100%</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">mr</span>
    
</div>';
		
		$calculations = '	<hr><h2>Dane:</h2>
		Masa substancji rozpuszczonej (ms): <b>'.$ms.'g</b>
						<br>Masa roztworu (mr): <b>'.$ms.'g + '.$input_array["mh2o"].'g = '.$mr.'</b>
						<br>Stężenie roztworu (Cp): <b>?</b>
												<hr><h2>Szukane:</h2>
						<br>Stężenie roztworu (Cp): <b>?</b>

						<hr><h2>Wzór:</h2>
						<br> <b>Cp = &nbsp;'.$frac1.'</b>
						<hr><h2>Obliczenia:</h2>
						<br>Cp = <div class="frac">
    <span>'.$ms.'g &bull;	 100%</span>
    <span class="symbol">/</span>
    <span class="bottom">'.$mr.'g</span>
</div> = <b>'.$score.'%</b>';
		
	}
	
	
	
	if ($input_array["template_id"] == "11") {
		
		///// USED VARIABLES ////////////////////////////////
		
		$weight_1 = $input_array["input_weight_1"];
		$weight_2 = $input_array["input_weight_2"];
		$element_1 = $input_array["input_element_1"];
		$element_2 = $input_array["input_element_2"];
		
		///////////////////////////////////////////////
		
		$element_1_parsed = $element_1;
		$element_2_parsed = $element_2;
		
		$atomic_mass_1 =  round (calculate_molar_mass ($element_1_parsed, $options_array = null),  $input_array["round_dec_num_atomic_mass"]);
		$atomic_mass_2 =  round (calculate_molar_mass ($element_2_parsed, $options_array = null),  $input_array["round_dec_num_atomic_mass"]);
		
		$atoms_element_1 = round ( ($weight_1 * MOL)/$atomic_mass_1,  $input_array["round_dec_num"]);
		$atoms_element_2 = round ( ($weight_2 * MOL)/$atomic_mass_2,  $input_array["round_dec_num"]);
		
		///////////////////////////////////////////////
		
		if ($atoms_element_1 > $atoms_element_2) {
			
		} elseif ($atoms_element_1 < $atoms_element_2) {
			
		} elseif ($atoms_element_1 == $atoms_element_2) {
			
		}
		
		///////////////////////////////////////////////
		
		$frac1 = '<div class="frac">
    <span>ms x 100%</span>
    <span class="symbol" style="background: white; color: white;">/</span>
    <span class="bottom">mr</span>
    
</div>';
		
		$calculations = '	<hr><h2>Dane:</h2>
		Masa substancji rozpuszczonej (ms): <b>'.$ms.'g</b>
						<br>Masa roztworu (mr): <b>'.$ms.'g + '.$input_array["mh2o"].'g = '.$mr.'</b>
						<br>Stężenie roztworu (Cp): <b>?</b>
												<hr><h2>Szukane:</h2>
						<br>Stężenie roztworu (Cp): <b>?</b>

						<hr><h2>Wzór:</h2>
						<br> <b>Cp = &nbsp;'.$frac1.'</b>
						<hr><h2>Obliczenia:</h2>
						<br>Cp = <div class="frac">
    <span>'.$ms.'g &bull;	 100%</span>
    <span class="symbol">/</span>
    <span class="bottom">'.$mr.'g</span>
</div> = <b>'.$score.'%</b>';
		
	}
	
	
	
	$output_array["solution_score"] = $score;
	$output_array["solution_calculations"] = $calculations;
	
	
	return $output_array;
	
}


function get_hydrocarbon_type ($hydrocarbon, $options_array) {
	
		$atoms_numbers_array = get_atom_number_in_compound ($hydrocarbon, $options_array);
	
		$hydrogen = $atoms_numbers_array["H"]; /// check number of hydrogen atoms
		$carbon = $atoms_numbers_array["C"]; /// check number of carbon atoms
		
		//// CHECK USING REGULAR FORMULA /////
		$for_alkan = 2*$carbon+2;
		$for_alken = 2*$carbon;
		$for_alkyn = 2*$carbon-2;
		
		if ($hydrocarbon == $for_alkan) {
			$hydrocarbon_type = "Alkan";
			
		} elseif ($hydrocarbon == $for_alken) {
			$hydrocarbon_type = "Alken";
			
		} elseif ($hydrocarbon == $for_alkyn) { /// alkin
			$hydrocarbon_type = "Alkin";
			
		} 
		
		return $hydrocarbon_type;
	
}

/// combustion reaction ///

function do_combustion_reaction ($input_array, $options_array)
{

	// // USING REGULAR FORMULA ///////
	// // ALKANES, ALKENES, ALKYNES ///

	if ($input_array["combusted_compound_type"] == "hydrocarbon") { /// sprawdzamy czy to jest węglowodór
		$combusted_compound = $input_array["combusted_compound"];
		$atoms_numbers_array = get_atom_number_in_compound($combusted_compound, $options_array);
		$n = $atoms_numbers_array["C"]; /// check number of carbon atoms
		$hydrocarbon_type = get_hydrocarbon_type($combusted_compound, $options_array);
		if ($input_array["combustion_type"] == "complete") { ///// spalanie całkowite
			if ($input_array["output_type"] == "detailed") {
				if ($hydrocarbon_type == "alkan") {
					$oxygen_compouding_number = (3 * $n + 1);
					$water_compouding_number = $n + 1;
					if ($oxygen_compouding_number % 2 == 0) {
						$oxygen_compouding_number = $oxygen_compouding_number / 2;
						$reaction_ready_array[0] = $combusted_compound . " + " . $oxygen_compouding_number . "O<sub>2</sub>" . " &rarr; " . $n . "CO<sub>2</sub>" . " + " . $water_compouding_number . "H<sub>2</sub>O";
					}
					else {
						$fraction1 = '<div class="frac"><span>' . $oxygen_compouding_number . '</span><span class="symbol" style="background: white; color: white;">/</span><span class="bottom">2</span></div>';
						$n2 = $n * 2;
						$water_compouding_number_double = $water_compouding_number * 2;
						$reaction_ready_array[0] = $combusted_compound . " + " . $fraction1 . "O<sub>2</sub>" . " &rarr; " . $n . "CO<sub>2</sub>" . " + " . $water_compouding_number . "H<sub>2</sub>O  &nbsp; / &bull;2";
						$reaction_ready_array[1] = "2" . $combusted_compound . " + " . $oxygen_compouding_number . "O<sub>2</sub>" . " &rarr; " . $n2 . "CO<sub>2</sub>" . " + " . $water_compouding_number_double . "H<sub>2</sub>O";
					}
				}
			}
			else {
				$unbalanced_reaction_chain = $combusted_compound . " + " . "O2 = CO2 + H2O";
				$reaction_balanced_array = balance_reaction($unbalanced_reaction_chain);
			}
		}

		if ($input_array["combustion_type"] == "semi_combustion") { ///// półspalanie
			if ($input_array["output_type"] == "detailed") {
			}
			else {
				$unbalanced_reaction_chain = $combusted_compound . " + " . "O2 = CO + H2O";

				//	echo $unbalanced_reaction_chain;

				$reaction_balanced_array = balance_reaction($unbalanced_reaction_chain);
			}
		}

		if ($input_array["combustion_type"] == "incomplete") { ///// spalanie całkowite
			if ($input_array["output_type"] == "detailed") {
			}
			else {
				$unbalanced_reaction_chain = $combusted_compound . " + " . "O2 = C + H2O";
				$reaction_balanced_array = balance_reaction($unbalanced_reaction_chain);
			}
		}
	}

	if ($input_array["output_type"] != "detailed") {
		$output_array["products"] =  $reaction_balanced_array[1];
		$output_array["substrats"] =  $reaction_balanced_array[0];
		return $output_array;
		//return $reaction_balanced_array[0] . " = " . $reaction_balanced_array[1];
	}
}


////////////////// DETERMINE COMPOUND STRUCTURAL FORMULA FROM POLISH //////////////////////////

function determine_compound_structural_formula ($chemical_compound_name) {
	
	
}

//////////////////////////// GET COMPOUND POLISH NAME FROM STRUCTURAL FORMULA //////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

function get_salt_force_array ($chemical_compound, $options_array) {
	
	global $forceAcidMoiety;
	
	$compound_type = get_compound_type($chemical_compound);
	
	$atoms_array = divide_to_atoms($chemical_compound , true, $options_array);
	
	if ($compound_type == "salt") {
	foreach ($atoms_array as $value) {
	
	
		if (get_compound_type($value) == "metal") {
			//// if checking alcali force
			$element_group = get_chemical_element_property ("group", $value);
			
			if ($element_group == 1 || $element_group == 2) {
				$element_force["alcali"] = "strong";
			} else {
				$element_force["alcali"] = "weak";
			}
			
		} elseif (contains_acid_moiety($value)) {
		   //// if checking acid force
			if ($forceAcidMoiety[$value] == 1) {
				$element_force["acid"] = "strong";
			} else {
				$element_force["acid"] = "weak";
			}
			
		}
		
	
	
	}
	
	} else {
		$output_array["errors"][] = "Compound isn't salt!";
	}
	
	$output_array["force_array"] = $element_force;
	
	
	return $output_array;
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////

function get_hydrolysis_reaction_type ($input_chemical_compound, $options_array) {
	
	$force_output = get_salt_force_array ($input_chemical_compound, $options_array);
	$force_array = $force_output["force_array"];
	
	if ($force_array["alcali"] == "weak" && $force_array["acid"] == "strong") {
		$type =  "Cationic";
		$odczyn =  "Basic";
	} elseif ($force_array["alcali"] == "strong" && $force_array["acid"] == "weak") {
		$type = "Anionic";
		$odczyn = "Acid";
	} elseif ($force_array["alcali"] == "weak" && $force_array["acid"] == "weak")  {
		$type = "Cationic-anionic";
		$odczyn = "Inert";
	} else {
		return false;
	}
	
	$output_array["recation"] = $odczyn;
	$output_array["type"] = $type;
	
	return $output_array;
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////

function translate_strings ($input, $lang) {
	global $TRANSLATION;
	
	if (is_array($input)) {
		foreach ($input as &$value) {
			$value = $TRANSLATION[$lang][$input];
		}
		return $value;
	} else {
		return $TRANSLATION[$lang][$input];
	}
	
}


//////////////////////////////////////////////////////////////////////////////////////////////////

function deobfuscate_chemical_compound ($input_chemical_compound, $options_array) {
	
	if ($input_chemical_compound == "ClH") {
		$output = "HCl";
	} else {
		$output = $input_chemical_compound;
	}
	
	
	return $output;
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////

function get_specified_atoms_number_in_compound ($chemical_compound, $specified_chemical_element, $options_array) {
	
	if (strpos($chemical_compound, $specified_chemical_element) !== false ) {
	
		$atoms_numbers_array = get_atom_number_in_compound ($chemical_compound, $options_array);
	
		$atoms_count = $atoms_numbers_array[$specified_chemical_element]; 
		
		return $atoms_count;
	
	} else {
		return false;
	}
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////

function do_hydrolysis_reaction($input_array, $options_array)
{

	// / REACTION SCHEME: AB + H2O ⇌ BH + AOH ////////
	// // SECTION FOR SALTS //////

	
	if ($input_array["hydrolysis_reaction_compound_type"] == "salt") {
		$force_output = get_salt_force_array($input_array["chemical_compound"], $options_array_force_check);
		$force_array = $force_output["force_array"];
		if ($force_array["alcali"] == "strong" && $force_array["acid"] == "strong") {
			$errors[] = "Reakcja hydrolizy nie zachodzi w przypadku soli mocnego kwasu i mocnej zadady!";
		}
		else {
			$atoms_array = divide_to_atoms($input_array["chemical_compound"], true, $options_array);
			$b = $atoms_array[0];
			$a = $atoms_array[1];
			if (is_non_oxysalt($input_array["chemical_compound"])) {
				if ($options_array["use_default_valence"] != 1) {
					$b_atoms = get_specified_atoms_number_in_compound($input_array["chemical_compound"], $b, $options_array);
					if ($b_atoms == 1) {
						$b_atoms = "";
					}

					$b_ready = $b . $b_atoms;
					$a_atoms = get_specified_atoms_number_in_compound($input_array["chemical_compound"], $a, $options_array);
					if ($a_atoms == 1) {
						$a_atoms = "";
					}

					$a_ready = $a . $a_atoms;
					$substrat_string = $b_ready . $a_ready;
					$product1_string = "H" . $b_atoms . $a;
					if ($a_atoms == 1 || $a_atoms == "") {
						$product2_string = $b . "OH";
					}
					else {
						$product2_string = $b . "(OH)" . $a_atoms;
					}
				}
				else {
					$substrat_string = deobfuscate_chemical_compound(order_compound_atoms_from_array(apply_valence($b . $a, true)));
					$product1_string = deobfuscate_chemical_compound(order_compound_atoms_from_array(apply_valence($a . "H", true)));
					$product2_string = deobfuscate_chemical_compound(order_compound_atoms_from_array(apply_valence($b . "OH", true)));
				}
			}
			else {
				if ($options_array["use_default_valence"] != 1) {
					$b_atoms = get_specified_atoms_number_in_compound($input_array["chemical_compound"], $b, $options_array);
					if ($b_atoms == 1) {
						$b_atoms = "";
						$b_ready = $b;
					}
					else {
						$b_ready = "(" . $b . ")" . $b_atoms;
					}

					$a_atoms = get_specified_atoms_number_in_compound($input_array["chemical_compound"], $a, $options_array);
					if ($a_atoms == 1) {
						$a_atoms = "";
					}

					$a_ready = $a . $a_atoms;
					$substrat_string = $a_ready . $b_ready;
					$product1_string = "H" . $a_atoms . $b;
					if ($b_atoms == 1 || $b_atoms == "") {
						$product2_string = $a . "OH";
					}
					else {
						$product2_string = $a . "(OH)" . $b_atoms;
					}
				}
				else {
					$substrat_string = deobfuscate_chemical_compound(order_compound_atoms_from_array(apply_valence($b . $a, true)));
					$product1_string = deobfuscate_chemical_compound(order_compound_atoms_from_array(apply_valence($b . "H", true)));
					$product2_string = deobfuscate_chemical_compound(order_compound_atoms_from_array(apply_valence($a . "OH", true)));
				}
			}

			// echo $b_ready;
			// echo $a;

			$internal_options_array["single_atoms_as_acid_moiety"] = 1; /// enable this option for valence appluing
			$unbalanced_reaction = $substrat_string . " + H2O = " . $product1_string . " + " . $product2_string;
		}
	}

	$balanced_array = balance_reaction($unbalanced_reaction);
	$ion_dissociation_options_array["in_water_valence"] = 1;
	$dissociated_array_0 = ion_dissociate($balanced_array[0], $type, $ion_dissociation_options_array);
	$dissociated_array_1 = ion_dissociate($balanced_array[1], $type, $ion_dissociation_options_array);
	echo $unbalanced_reaction;
	print_r($dissociated_array_1);
	echo $balanced_array[0] . " ⇌ " . $balanced_array[1];
	print_r($errors);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function calculate_molar_mass ($chemical_compound, $options_array = null) {
	
	if ($options_array["algorithm"] == 1) {
		
		$chemical_compound = remove_compounding_number($chemical_compound);
		$pathToRhinoJar = '/var/www/html/reactionsolver/engine/rhino1.7.7.1/lib/rhino-1.7.7.1.jar';
		$javascriptFile = '/var/www/html/reactionsolver/engine/stechiometria1.js';
		$output = shell_exec("java -jar $pathToRhinoJar $javascriptFile ".'"'.$chemical_compound.'"');
 
		return $output;
		
	} elseif ($options_array["algorithm"] == 2) {
		
		$chemical_compound = remove_compounding_number($chemical_compound);
		$javascriptFile = '/var/www/html/reactionsolver/engine/stechiometria1_phantom.js';
		$output = shell_exec("QT_QPA_PLATFORM=offscreen phantomjs $javascriptFile ".'"'.$chemical_compound.'"');
 
		return $output;
		
	}

}


///////////////////////////////////////////////////////////

function get_molar_mass_percent_array ($chemical_compound, $sum_molar_mass = null, $options_array = null) {
	
	$options_array["divide_atoms_include_special_compounds"] = true;
	$internal_options_array["algorithm"] = 2;
	$atoms_array = divide_to_atoms($chemical_compound , false, $options_array);
	
	if ($sum_molar_mass == null) {
		
		$sum_molar_mass = calculate_molar_mass ($chemical_compound, $internal_options_array);
		
	} 
	
	foreach ($atoms_array as $atom) {
			$multiliper = get_specified_atoms_number_in_compound ($chemical_compound, $atom, $options_array);
			$atoms_molar_mass_array_tmp[$atom] = calculate_molar_mass ($atom, $internal_options_array)*$multiliper;
			$atoms_molar_mass_array[$atom] = $atoms_molar_mass_array_tmp[$atom]/$sum_molar_mass;
			$atoms_molar_mass_percent_array[$atom] = $atoms_molar_mass_array[$atom]*100;
			//echo $multiliper;
	}
	
	return $atoms_molar_mass_percent_array;
	//return $sum_molar_mass;
	//return $atoms_arry;
	
	
}

///////////////////////////////////////////////////////////

function get_hydrocarbon_name ($hydrocarbon_formula, $options_array) {
	
	
	
	$atoms_numbers_array = get_atom_number_in_compound ($hydrocarbon_formula, $options_array);
	
	$hydrogen = $atoms_numbers_array["H"]; /// check number of hydrogen atoms
	$carbon = $atoms_numbers_array["C"]; /// check number of carbon atoms
	
	$hydrocarbon_type = get_hydrocarbon_type ($hydrocarbon_formula, $options_array);
	
	if ($options_array["lang"] = "PL") {
	
	if ($carbon == 1) {
		$output_name_tmp_1 = "met";
	} elseif ($carbon == 2) {
		$output_name_tmp_1 = "et";
	}  elseif ($carbon == 3) {
		$output_name_tmp_1 = "prop";
	}  elseif ($carbon == 4) {
		$output_name_tmp_1 = "but";
	}  elseif ($carbon == 5) {
		$output_name_tmp_1 = "pent";
	}  elseif ($carbon == 6) {
		$output_name_tmp_1 = "heks";
	}  elseif ($carbon == 7) {
		$output_name_tmp_1 = "hept";
	}  elseif ($carbon == 8) {
		$output_name_tmp_1 = "okt";
	}  elseif ($carbon == 9) {
		$output_name_tmp_1 = "non";
	}  elseif ($carbon == 10) {
		$output_name_tmp_1 = "dek";
	}
	
	if ($hydrocarbon_type == "Alkan") {
		$output_name_tmp_2 = $output_name_tmp_1."an" ;
	} elseif ($hydrocarbon_type == "Alken") {
		$output_name_tmp_2 = $output_name_tmp_1."en" ;
	} elseif ($hydrocarbon_type == "Alkin") {
		$output_name_tmp_2 = $output_name_tmp_1."yn" ;
	}
	
	}
	
}

///////////////////////////////////////////////////////////

function trim_all($str) {
return preg_replace('/\s/', "", $str);
}

/////////////////////////////////////////////////////////




?>
