data_dir(label::String) :: String = "raw_data/" * label * "/"
pp_dir(label::String) :: String = "pretty_things/" * label * "/"

function check_pp_dir(label::String)
    if !isdir(pp_dir(label))
        mkdir(pp_dir(label))
    end
end