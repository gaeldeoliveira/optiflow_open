classdef lg
    %LG is a simple static class collecting methods used for logging 
    % system messages from the Optiflow airfoil optimization system and the
    % VI calibration framework.
    %
    % File output is implemented in a simple way!
    %
    
    properties
        
    end
    
    methods
        function logger_object = lg()
            %Dummy constructor
        end
    end
    
    methods(Static)
        function msg(obj, message)
           % Ouput message to standard output and log file
           %  obj is the handle to the object/class emitting the message
           %  message is a string containing the message that should be 
           
           % Make message string
           msg_str = ['STD:' , datestr(now, 30) , ':' , class(obj) , '> ' , message];
           % Display message string
           disp(msg_str);
           % Write   message string to file
           fid = fopen('lg.log','a'); fprintf(fid,[msg_str,'\n']); fclose(fid);
        end
        
        function err(obj, message)
            % Output message to standard output and error file
            msg_str = ['ERR:' , datestr(now, 30) , ':' , class(obj) , '> ' , message];
            disp(msg_str);
            fid = fopen('err.log','a'); fprintf(fid,[msg_str,'\n']); fclose(fid);
        end
        
        function clr(message)
            % Move current log files   (ok, this is weird coding... like, really silly)
            flag_lg_dir_missing  = system(['cp lg.log ./tmp/lg/lg_upto_' datestr(now(), 30), '.log']);
            if flag_lg_dir_missing
                system('mkdir ./tmp/lg');
                system(['cp lg.log ./tmp/lg/lg_upto_' datestr(now(), 30), '.log']);
            end
            % Move current error files (ok, this is weird coding... like, really silly)
            flag_err_dir_missing = system(['cp err.log ./tmp/err/err_upto_' datestr(now(), 30), '.log']);
            if flag_err_dir_missing
                system('mkdir ./tmp/err');
                system(['cp err.log ./tmp/err/err_upto_' datestr(now(), 30), '.log'])
            end
            % Remove current files 
            system('rm -f err.log lg.log');
            % Start anew and say it out loud!
            lg.msg(lg,  '---- START NEW LOG BOOK AFTER REQUEST ----');
            lg.msg(lg, ['---- request issuer: ', message]);
            lg.err(lg,  '---- START NEW ERR BOOK AFTER REQUEST ----');
            lg.err(lg, ['---- request issuer: ', message]);
        end
    end
end

