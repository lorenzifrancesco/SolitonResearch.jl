
function colorof(key::Tuple{String, Float64})
  name = key[1]
  return colorof(name)
end

function colorof(name::String)
  if name == "G3"
    return :red
  elseif name == "G1"
    return :grey
  elseif name == "N"
    return :green 
  elseif name == "Np"
    return :green
  elseif name == "CQ"
    return :blue
  else 
    return :black
  end
end

function lineof(key::Tuple{String, Float64})
  name = key[1]
  return lineof(name)
end


function lineof(name::String)
  if name == "G3"
    return :solid
  elseif name == "G1"
    return :solid
  elseif name == "N"
    return :dot
  elseif name == "Np"
    return :dash
  elseif name == "CQ"
    return :dashdot
  else 
    return :solid
  end
end

function nameof(key::Tuple{String, Float64})
  name = key[1]
  return nameof(name)
end

function nameof(name::String)
  if name == "G3"
    return "3D-GPE"
  elseif name == "G1"
    return "1D-GPE"
  elseif name == "N"
    return "NPSE"
  elseif name == "Np"
    return "NPSE+"
  elseif name == "CQ"
    return "CQ-GPE"
  else 
    return "Unknown"
  end
end