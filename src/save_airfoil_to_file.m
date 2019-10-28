function [written_coordinates] = save_airfoil_to_file(SW, x, output_file)
    % Get system context
    SC = SW.SC;
    % Clean folder to write airfoil coordinates
    warning('Off', 'MATLAB:DELETE:FileNotFound')
    SW.clean_folders(1);
    warning('On' , 'MATLAB:DELETE:FileNotFound')
    % Compute and write airfoil polar coordinates
    %fval = GM.call_cost_function(x);
    coordinates_block_cell = SW.write_airfoils_block({x});
    % Extract written coordinates from cell
    written_coordinates = coordinates_block_cell{1};
    % Write path of airfoil that was just written
    airfoil_source_file = [SC.tmp_subdir, SC.core_subdirs{1} , SC.airfoil_filename];
    % Build string to copy airfoil filename
    airfoil_string_copy = ['cp ' , airfoil_source_file , ' ' ,output_file]; 
    % Copy airfoil to destination 
    f = system(airfoil_string_copy);
    % Inform of writting
    if (f==0)
        disp(['Wrote airfoil to: ' , output_file])
    end
    % Return!
end
