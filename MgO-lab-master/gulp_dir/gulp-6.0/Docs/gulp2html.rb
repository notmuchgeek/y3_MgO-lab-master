#!/usr/bin/ruby -w
# copyright Andrew Rohl 2015
# To run:
# ./gulp2html.rb help.txt > help.html

# TODO implement special formatting for space groups
 
help_filename = ARGV[0]
help_file = File.new(help_filename, "r")
keyword_hash = Hash.new{|hash, key| hash[key] = Array.new;}
option_hash = Hash.new{|hash, key| hash[key] = Array.new;}
info_hash = Hash.new{|hash, key| hash[key] = Array.new;}
lineno = 0
while line = help_file.gets
  lineno += 1
  line.chomp!
  if line.start_with?("@@")
    keyword = line[2..-1]
    # chomp next line as just keyword again
    help_file.gets
    lineno += 1
    # get type
    if keyword == "valid_spacegroups"
      type = "Information"
    else
      line = help_file.gets
      lineno += 1
      rawtype = line.split(":")[1].strip!
      rawtype.capitalize!
      type = rawtype.split(" ")[0]
      if type == "Keyword"
        keyword_hash[keyword].push line
      elsif type == "Option"
        option_hash[keyword].push line
      elsif type == "Information"
        info_hash[keyword].push line
      else
        puts "unknown type %s on line %d" % [type, lineno]
      end
    end
  else
    # skip blank lines
    if !line.empty?
      if type == "Keyword"
        keyword_hash[keyword].push line
      elsif type == "Option"
        option_hash[keyword].push line
      elsif type == "Information"
        info_hash[keyword].push line
      end
    end
  end
end

puts "<!DOCTYPE html>"
puts "<html>"
puts "<head>"
puts "<title>GULP help file</title>"
puts "<h1>GULP help file</h1>"
puts "</head>"

puts "<body>"
puts '<hr>'
puts '<h2>Keywords in alphabetical order</h2>'

help_hash = keyword_hash
help_hash.sort.map do|keywrd, help_text|
  puts "<a href=\#%s>%s</a>&nbsp;" % [keywrd, keywrd]
end

puts '<hr>'
puts '<h2>Options in alphabetical order</h2>'

help_hash = option_hash
help_hash.sort.map do|keywrd, help_text|
  puts "<a href=\#%s>%s</a>&nbsp;" % [keywrd, keywrd]
end

puts '<hr>'
puts '<h2>Information in alphabetical order</h2>'

help_hash = info_hash
help_hash.sort.map do|keywrd, help_text|
  puts "<a href=\#%s>%s</a>&nbsp;" % [keywrd, keywrd]
end
puts '<hr>'

hashes = [keyword_hash, option_hash, info_hash]
hashes.each do|hlp_hash| 
	hlp_hash.sort.map do|keywrd, help_text|
		puts "<a name=\"%s\">" % keywrd
		puts "<h3>%s</h3>" % keywrd
		puts "<table>"
		if (keywrd == "valid_spacegroups") || (keywrd == "trajectory_format")
			# special code to deal with space groups list/trajectory format
			puts "<pre>"
			help_text.each do |lne|
				puts lne
			end
			puts "</pre>"
		else
			have_see_also = false
			ignore_next_line = false
			help_text.each do |lne|
				if (ignore_next_line)
					ignore_next_line = false
				else
					puts "<tr>"
					# escape the brackets
					lne.gsub! "<", "&lt;"
					lne.gsub! ">", "&gt;"
					# each line is actually two columns - split on the first colon
					# hack - replace first colon with > and then split on it
					lne.sub! ":", ">"
					columns = lne.split ">"
					col0 = columns[0]
					# testing to see if we have hit the See also section of the keyword as we have to
					# add the links. Note that assume this is last part of keyword help
					if col0 == "See also"
						have_see_also = true
					elsif col0 == "and"
						have_see_also = false
						# ignore next line as just keyword again
						ignore_next_line = true
					end
					if col0.nil?
						puts "  <td></td>"
					elsif col0.strip.empty?
						puts "  <td></td>"
					else
						puts "  <td style=\"padding-right: 50px;\"><strong>" + col0 + "</strong></td>"
					end
					col1 = columns[1]
					if col1.nil?
						puts "  <td></td>"
					else
						col1.strip!
						if have_see_also
							puts "  <td>"
							keyword_to_ref = col1.split ","
							keyword_to_ref.each do |word|
								# remove any periods in the see also list
								word.sub! ".", ""
								puts "<a href=\#%s>%s</a>" % [word.strip, word.strip]
							end
							puts "</td>"
						else  
							puts "  <td>" + col1 + "</td>"
						end
					end
					puts "</tr>"
				end
			end
			puts "</table>"
			puts '<hr>'
		end
	end
end
puts "</body>"
puts "</html>"
