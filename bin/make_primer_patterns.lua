#!/usr/bin/env lua

-- Usage:
-- lua make_primer_patterns.lua -i input.fasta [-o output_prefix] [-f forward_pattern] [-r reverse_pattern]

local function parse_args()
	local opts = {
		input_fasta = nil,
		output_prefix = "primer_patterns",
		forward_pattern = "^(.*?)",
		reverse_pattern = "^(.*?)",
	}

	local i = 1
	while i <= #arg do
		local a = arg[i]
		if a == "-i" or a == "--input_fasta" then
			i = i + 1
			opts.input_fasta = arg[i]
		elseif a == "-o" or a == "--output_prefix" then
			i = i + 1
			opts.output_prefix = arg[i]
		elseif a == "-f" or a == "--forward_pattern" then
			i = i + 1
			opts.forward_pattern = arg[i]
		elseif a == "-r" or a == "--reverse_pattern" then
			i = i + 1
			opts.reverse_pattern = arg[i]
		else
			error("Unknown argument: " .. a)
		end
		i = i + 1
	end

	assert(opts.input_fasta, "You must provide --input_fasta (-i)")
	return opts
end

local function read_fasta_lines(path)
	local file, err = io.open(path, "r")
	assert(file, "Could not open FASTA file: " .. (err or "unknown error"))
	local lines = {}
	for line in file:lines() do
		lines[#lines + 1] = line:gsub("%s+", "") -- strip whitespace
	end
	file:close()
	return lines
end

local function generate_patterns(fasta_path, label, fwd_prefix, rev_suffix)
	local lines = read_fasta_lines(fasta_path)

	-- Extract sequences and headers
	local seqs, headers, entries = {}, {}, {}
	local accumulator = nil
	for _, line in ipairs(lines) do
		if line:sub(1, 1) == ">" then
			if accumulator then
				entries[#entries + 1] = accumulator
			end
			accumulator = { header = line, seq = "" }
		else
			assert(accumulator)
			accumulator.seq = accumulator.seq .. line
		end
	end
	if accumulator then
		entries[#entries + 1] = accumulator
	end

	-- Crash if the number of parsed sequences isn't exactly 2
	assert(#seqs == 2)

	-- Parse start coordinates from header lines
	local starts = {}
	for _, header in ipairs(headers) do
		local start = header:match(":(%d+)%-%d+")
		assert(start, "Invalid header format (expected bedtools-style): " .. header)
		starts[#starts + 1] = tonumber(start)
	end

	-- Heuristic check for orientation assumption
	if starts[1] > starts[2] then
		io.stderr:write(
			"⚠️  Warning: Please double check that the provided FASTA is formatted like an output from `bedtools getfasta`, e.g.\n\n'>PP599462.1:0-16'"
		)
	end

	-- Build patterns
	local fwd_pattern = fwd_prefix .. seqs[1]
	local rev_pattern = seqs[2] .. rev_suffix

	-- Write output
	local out, err = io.open(label .. ".txt", "w")
	assert(out, "Failed to write output: " .. (err or "unknown error"))
	out:write(fwd_pattern .. "\n")
	out:write(rev_pattern .. "\n")
end

local function main()
	local opts = parse_args()

	-- make sure the provided file exists
	local file = io.open(opts.input_fasta, "r")
	assert(file, "Input FASTA file does not exist: " .. opts.input_fasta)
	file:close()

	-- generate the patterns and write them to a text file
	generate_patterns(opts.input_fasta, opts.output_prefix, opts.forward_pattern, opts.reverse_pattern)
end

-- Run main
if debug.getinfo(1, "S").short_src == arg[0] then
	main()
end
