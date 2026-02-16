using Dates

function doubleoutput(output::String, label::String)
    output = "[" * string(now()) * "]  " * output * "\n"
    outputfilename = "raw_data/" * label * "/log.txt"
    file = open(outputfilename, "a")
    write(file, output)
    close(file)
    print(output)
end