#! /usr/bin/env ruby

#version 2, 2009/Jan/24
#version 3, 2015/Dec/8


if ARGV.length == 1
	scale = 1.0
elsif ARGV.length == 2
	scale = ARGV.shift.to_f
else
	STDERR.puts "USAGE: newick2mafft.rb scale input_tree > output"
	exit
end

if scale <= 0.0 then
	STDERR.puts "Inappropriate scale, #{scale.to_s}"
	exit
end

STDERR.puts "scale = " + scale.to_s

infp = File.open( ARGV.shift, "r" )

tree = ""
while line = infp.gets
	tree += line.strip
	break if tree =~ /;$/
end
infp.close


#tree = tree.gsub( /_.*?:/, ":" ).gsub(/[0-9]\.[0-9]*e-[0-9][0-9]/, "0").gsub(/\[.*?\]/,"").gsub(/ /, "").gsub(/:\-[0-9\.]+/, ":0.0" )
#tree = tree.gsub( /_.*?:/, ":" ).gsub(/[0-9]\.[0-9]*e-[0-9][0-9]/, "0").gsub(/\[.*?\]/,"").gsub(/ /, "")
tree = tree.gsub( /_.*?:/, ":" ).gsub(/[0-9\.]*[eE]-[0-9]*/, "0").gsub(/\[.*?\]/,"").gsub(/ /, "")


STDERR.puts "Initial tree = " +  tree

def resolve( tree )


while 1
#	p tree
	tree.sub!( /\,([0-9]+):(\-?[0-9\.]+)\,([0-9]+):(\-?[0-9\.]+)/, ",XXX" )
	hit1 = $1
	hit2 = $2
	hit3 = $3
	hit4 = $4

#	p hit1
#	p hit2
#	p hit3
#	p hit4

#	puts "introduce XXX"
#	p tree

	break unless tree.index(/XXX/)

	poshit = tree.index(/XXX/)
#	puts "poshit=" + poshit.to_s

	i = poshit
	height = 0
	while i >= 0
		break if height == 0 && tree[i..i] == '('
		if tree[i..i] == ')' then
			height += 1
		elsif tree[i..i] == '(' then
			height -= 1
		end
		i -= 1
	end

	poskakko = i
#	puts "poskakko = " + poskakko.to_s
	zenhan = tree[0..poskakko]
	zenhan = "" if poskakko == -1
#	puts "zenhan = " + zenhan

	treelen = tree.length
	tree = zenhan + "(" + tree[poskakko+1..treelen]
#	puts "add ("
#	p tree
	tree.sub!( /XXX/, "#{hit1}:#{hit2}):0,#{hit3}:#{hit4}" )

#	p tree
end


return tree

end

memi = [-1,-1]
leni = [-1,-1]

while tree.index( /\(/ ) 

	tree = resolve( tree )

	tree.sub!( /\(([0-9]+):(\-?[0-9\.]+),([0-9]+):(\-?[0-9\.]+)\)/, "XXX" )
	memi[0] = $1.to_i
	leni[0] = $2.to_f * scale
	memi[1] = $3.to_i
	leni[1] = $4.to_f * scale

	if leni[0] > 10 || leni[1] > 10 then
		STDERR.puts ""
		STDERR.puts "Please check the scale of branch length!"
		STDERR.puts "The unit of branch lengths must be 'substitution/site'"
		STDERR.puts "If the unit is 'substition' in your tree, please"
		STDERR.puts "use the scale argument,"
		STDERR.puts "% newick2mafft scale in > out"
		STDERR.puts "where scale = 1/(alignment length)"
		STDERR.puts ""
		exit 1
	end

#	STDERR.puts "subtree = " + $&

	if memi[1] < memi[0] then
		memi.reverse!
		leni.reverse!
	end

	tree.sub!( /XXX/, memi[0].to_s )

#	STDERR.puts "Tree = " + tree

	printf( "%5d %5d %10.5f %10.5f\n", memi[0], memi[1], leni[0], leni[1] )

end
