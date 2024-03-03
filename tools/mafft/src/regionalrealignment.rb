#! /usr/bin/env ruby

$MAFFTCOMMAND = '"/usr/local/bin/mafft"'
# Edit the above line to specify the location of mafft.
# $MAFFTCOMMAND = '"C:\folder name\mafft.bat"' # windows
# $MAFFTCOMMAND = '"/usr/local/bin/mafft"'     # mac or cygwin
# $MAFFTCOMMAND = '"/usr/bin/mafft"'           # linux (rpm)
# $MAFFTCOMMAND = '"/somewhere/mafft.bat"'     # all-in-one version for linux or mac

#####################################################################
#
# regionalrealignment.rb version 0.2 (2013/Sep/21)
# ruby regionalrealignment.rb setting input > output
# See http://mafft.cbrc.jp/alignment/software/regionalrealignment.html
# 
# v0.2, 2013/Sep/21, Fixed a windows-specific bug.
#
#####################################################################


def readfasta( fp, name, seq )
        nseq = 0
        tmpseq = ""
        while fp.gets
                if $_ =~ /^>/ then
                        name.push( $_.sub(/>/,"").strip )
                        seq.push( tmpseq ) if nseq > 0
                        nseq += 1
                        tmpseq = ""
                else
                        tmpseq += $_.strip
                end
        end
        seq.push( tmpseq )
        return nseq
end

def resolve( tree )
	while 1
#		p tree
		tree.sub!( /\,([0-9]+):(\-?[0-9\.]+)\,([0-9]+):(\-?[0-9\.]+)/, ",XXX" )
		hit1 = $1
		hit2 = $2
		hit3 = $3
		hit4 = $4
	
#		p hit1
#		p hit2
#		p hit3
#		p hit4
	
#		puts "introduce XXX"
#		p tree
	
		break unless tree.index(/XXX/)
	
		poshit = tree.index(/XXX/)
#		puts "poshit=" + poshit.to_s
	
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
#		puts "poskakko = " + poskakko.to_s
		zenhan = tree[0..poskakko]
		zenhan = "" if poskakko == -1
#		puts "zenhan = " + zenhan
	
		treelen = tree.length
		tree = zenhan + "(" + tree[poskakko+1..treelen]
#		puts "add ("
#		p tree
		tree.sub!( /XXX/, "#{hit1}:#{hit2}):0,#{hit3}:#{hit4}" )
	
#		p tree
end


return tree

end

if ARGV.length != 2 then
	STDERR.puts ""
	STDERR.puts "Usage: ruby #{$0} setingfile inputfile > output"
	STDERR.puts ""
	exit 1
end

infilename = ARGV[1]
tname = []
tseq = []
infp = File.open( infilename, "r" )
tin = readfasta( infp, tname, tseq )
infp.close

if tin == 0 then
		STDERR.puts ""
		STDERR.puts "Error in the '#{infilename}' file.  Is this FASTA format?\n"
		STDERR.puts ""
		exit 1
end

alnlen = tseq[0].length
if alnlen == 0 then
		STDERR.puts ""
		STDERR.puts "Error in the '#{infilename}' file.  Is this FASTA format?\n"
		STDERR.puts ""
		exit 1
end


for i in 0..(tin-1)
	if alnlen != tseq[i].length then
		STDERR.puts ""
		STDERR.puts "Please insert gaps such that all the input sequences have the same length.\n"
		STDERR.puts ""
		exit 1
	end
end

checkmap = []
for i in 0..(alnlen-1)
	checkmap.push(0)
end

outputseq = []
for i in 0..(tin-1)
	outputseq.push("")
end


settingfile = ARGV[0].to_s
reg = []
startpos = []
endpos = []
realign = []
options = []
treeoption = ""
revwarn = 0
sfp = File.open( settingfile, "r" )
while line = sfp.gets
	line.sub!(/#.*/,"")
	next if line.length < 2
	if line.strip =~ /^treeoption / then
		treeoption = line.strip.sub(/.*treeoption/,"")
		break
	end
end
sfp.close
sfp = File.open( settingfile, "r" )
while line = sfp.gets
	line.sub!(/#.*/,"")
	next if line.length < 2
	next if line.strip =~ /^treeoption/
	startposv = line.split(' ')[0].to_i - 1
	endposv = line.split(' ')[1].to_i - 1
	if startposv < 0 || endposv < 0 then
		STDERR.puts "\nError in the '#{settingfile}' file. Please check this line:\n"
		STDERR.puts line
		STDERR.puts "Sites must be numbered as 1, 2, ...\n"
		STDERR.puts "\n"
		exit 1
	end
	if startposv >= alnlen || endposv >= alnlen then
		STDERR.puts "\nError in the '#{settingfile}' file. Please check this line:\n"
		STDERR.puts line
		STDERR.puts "Sites must be numbered as 1, 2, ... #{alnlen}\n"
		STDERR.puts "\n"
		exit 1
	end
	if startposv > endposv then
		STDERR.puts "\nWarning. Please check this line:\n"
		STDERR.puts line
		STDERR.puts "Start position > End position ?\n"
		STDERR.puts "\n"
		revwarn = 1
#		exit 1
	end
	startpos.push( startposv )
	endpos.push( endposv )
	if startposv > endposv
		for k in (endposv)..(startposv)
			checkmap[k] += 1
		end
	else
		for k in (startposv)..(endposv)
			checkmap[k] += 1
		end
	end
	if line.split(' ')[2] == "realign"  then
		realign.push( 1 )
	elsif line.split(' ')[2] == "preserve" then
		realign.push( 0 )
	else
		STDERR.puts "\n"
		STDERR.puts "The third column must be 'realign' or 'preserve'\n"
		STDERR.puts "Please check this line:\n"
		STDERR.puts line
		STDERR.puts "\n"
		exit 1
	end
	if line =~ / \-\-/ && line =~ /realign/ then
		options.push( line.sub(/.*realign/,"").strip )
	else
		options.push( treeoption )
	end
end
sfp.close

#p startpos
#p endpos
#p options


#res = system "#{$MAFFTCOMMAND} #{treeoption} --treeout --retree 0 --thread -1 #{infilename} > _dum"
res = system "#{$MAFFTCOMMAND} #{treeoption} --treeout --retree 0  #{infilename} > _dum"

if res == false then
	STDERR.puts "\n"
	STDERR.puts "ERROR in building a guide tree"
	STDERR.puts "\n"
	exit 1
end

treefp = File.open( "#{infilename}.tree", "r" )

tree = ""
while line = treefp.gets
	tree += line.strip
	break if tree =~ /;$/
end
treefp.close

tree = tree.gsub( /_.*?:/, ":" ).gsub(/[0-9]\.[0-9]*e-[0-9][0-9]/, "0").gsub(/\[.*?\]/,"").gsub(/ /, "")
scale = 1.0
mtreefp = File.open("_tree", "w")


#STDERR.puts "Tree = " +  tree

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

	mtreefp.printf( "%5d %5d %10.5f %10.5f\n", memi[0], memi[1], leni[0], leni[1] )

end


mtreefp.close

numreg = startpos.length

for i in 0..(numreg-1)

	partfp = File.open( "_part", "w" )
	for j in 0..(tin-1)
		partfp.puts ">" + tname[j]
		if startpos[i] > endpos[i] then
			partfp.puts tseq[j][endpos[i]..startpos[i]].reverse
		else
			partfp.puts tseq[j][startpos[i]..endpos[i]]
		end
	end
	partfp.close

	if( realign[i] == 1 ) then
		STDERR.puts "Aligning region #{startpos[i]+1} - #{endpos[i]+1}"
		res = system "#{$MAFFTCOMMAND} #{options[i]} --inputorder --treein _tree _part > _partout"
		if res == false then
			STDERR.puts "\n"
			STDERR.puts "ERROR in aligning region #{startpos[i]+1} - #{endpos[i]+1}"
			STDERR.puts "Please check the option:"
			STDERR.puts "#{options[i]}"
			STDERR.puts "\n"
			exit 1
		end

	else
		STDERR.puts "Copying region #{startpos[i]+1} - #{endpos[i]+1}"
#		system "cp _part _partout"
		File.rename( "_part", "_partout" )
	end

	pname = []
	pseq = []
	partfp = File.open( "_partout", "r" )
	pin = readfasta( partfp, pname, pseq )
	partfp.close
	for j in 0..(tin-1)
		outputseq[j] += pseq[j]
	end
end

for j in 0..(tin-1)
	puts ">" + tname[j]
	puts outputseq[j]
end

STDERR.puts "Done."

numdupsites = checkmap.select{|x| x>1}.length
if numdupsites > 0 then
	STDERR.puts ""
	STDERR.puts "#########################################################"
	STDERR.puts "# Warning: #{numdupsites} sites were duplicatedly selected."
	STDERR.puts "#########################################################"
	STDERR.puts ""
end

numunselectedsites = checkmap.select{|x| x==0}.length
if numunselectedsites > 0 then
	STDERR.puts ""
	STDERR.puts "#########################################################"
	STDERR.puts "# Warning: #{numunselectedsites} sites were not selected."
	STDERR.puts "#########################################################"
	STDERR.puts ""
end

if revwarn == 1 then
	STDERR.puts ""
	STDERR.puts "#########################################################"
	STDERR.puts "# Warning: The order of sites were reversed."
	STDERR.puts "#########################################################"
	STDERR.puts ""
end


STDERR.puts ""
STDERR.puts "           Tree: computed  with #{treeoption} --treeout "
for i in 0..(numreg-1)
	range = sprintf( "%6d - %6d", startpos[i]+1, endpos[i]+1 )
	if realign[i] == 1 then
		STDERR.puts "#{range}: realigned with #{options[i]} --treein (tree)"
	else
		STDERR.puts "#{range}: preserved"
	end
end
STDERR.puts ""

File.delete( "_dum" )
File.delete( "_tree" )
File.delete( "_part" )
File.delete( "_partout" )

