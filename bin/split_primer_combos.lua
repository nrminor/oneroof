#!/usr/bin/env lua

-- Usage: lua split_primers.lua input.bed _LEFT _RIGHT
-- Default suffixes if not provided
local input_bed = arg[1] -- remember lua is 1-based!
local forward_suffix = arg[2] or "_LEFT"
local reverse_suffix = arg[3] or "_RIGHT"

-- Make sure the input bed file is provided
assert(input_bed, "Usage: lua split_primers.lua <input.bed> [_LEFT] [_RIGHT]")

-- Initialize a table to store each primer
local primers = {}

-- Read BED line by line
local file = io.open(input_bed, "r")
assert(file, "Failed to open input BED file.")

for line in file:lines() do
	local ref, start_pos, stop_pos, name, index, sense =
		line:match("([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)")
	assert(name, "Malformed BED line: " .. line)

	local base_name = name
	base_name = base_name:gsub(forward_suffix, "")
	base_name = base_name:gsub(reverse_suffix, "")

	primers[base_name] = primers[base_name] or {}
	table.insert(primers[base_name], line)
end

file:close()

-- Write each group to its own BED file
for base_name, records in pairs(primers) do
	assert(#records == 2, "Expected 2 records for " .. base_name .. ", found " .. #records)
	local out = assert(io.open(base_name .. ".bed", "w"))
	for _, rec in ipairs(records) do
		out:write(rec .. "\n")
	end
	out:close()
end
