use strict;
use warnings;
use v5.10.0;
use FindBin;
use Data::Dumper;
use GD;
use POSIX;
use Storable qw(dclone);
use List::Util qw(min max sum);
use lib "$FindBin::Bin/src/";
use PathLayDBs qw(metaDB);
$Data::Dumper::Indent = 1;

# general settings
my $debug = 1;
my $organism_code = "hsa";
my $db = "kegg";
my $gmt_file = "${organism_code}.${db}.meta.gmt";
my $nodes_folder = "$FindBin::Bin/pathlay_data/${organism_code}/maps/${db}/";



# data loading for mapping
# metadata and cross ref
my $metadata_file = "metadata_lipids.tsv";
my $lipids = load_metadata(file => $metadata_file);

# de data (based on metadata content)
my $de_file = "de_lipids.tsv";
$lipids = load_de_data(file => $de_file, metadata => $lipids);

# kegg meta gmt loading
my $gmt_path = "$FindBin::Bin/pathlay_data/${organism_code}/db/${db}/${gmt_file}";
my $gmt = load_gmt(file => $gmt_path, metadata => $lipids);


my $pathways = {};
foreach my $pID (sort keys %{$gmt->{pathway2meta}}) {
	say "${pID} - $gmt->{metadata}->{id2name}->{$pID}"  if $debug;
	my $node_file = "${nodes_folder}${pID}.nodes";
	my @nodes = load_nodes(file => $node_file, data => $lipids);
	foreach my $node (@nodes) {
		load_de_on_node(node => $node, data => $lipids);
	}
	my $complexes_index = make_complexes(nodes => \@nodes);
	# say Dumper $complexes_index;


	# plot on imgs
	my $imgPath = "pathlay_data/pathways/${db}/";
	my $pImg = $pID;
	$pImg =~ s/${organism_code}//;
	open my $in, '<:raw', "${imgPath}${pImg}.png" or die "Cannot open ${pImg}: $!";
	my $bg = GD::Image->newFromPng($in) or die "Invalid PNG image";
	close $in;

	

	foreach my $coord (keys %$complexes_index) {
		my $counts = {
			dev_down => 0,
			dev_up => 0
		};
		foreach my $n (keys %{$complexes_index->{$coord}}) {
			foreach my $lipid (keys %{$complexes_index->{$coord}->{$n}->{lipids}}) {
				my $lipidObj = $complexes_index->{$coord}->{$n}->{lipids}->{$lipid};
				$counts->{dev_up}++ if ($lipidObj->{dev} > 0);
				$counts->{dev_down}++ if ($lipidObj->{dev} < 0);
			}
		}
		my @counts_dev = ($counts->{dev_up},$counts->{dev_down});
		my ($x,$y) = split(",",$coord);
		say("Complex Coords: x:$x - y:$y") if ($debug);

		my $graph_size = 100;
		my $graph = new GD::Image($graph_size,$graph_size);

		my $white = $graph->colorAllocate(255,255,255);
		my $red = $graph->colorAllocate(255, 0, 0);
		my $green = $graph->colorAllocate(0,255,0);
		my @colors_dev = ($red,$green);

		$graph->transparent($white);
		$graph->interlaced('true');
		$graph->filledRectangle(0, 0, $graph_size, $graph_size, $white);

		$graph = sliced_crown_arc(
			graph => $graph,
			centerX => $graph_size/2,
			centerY => $graph_size/2,
			diameter1 => 20,
			diameter2 => 30,
			counts => \@counts_dev,
			colors => \@colors_dev
		);

		# Add text (optional)
		# $bg->string(gdSmallFont, 60, 110, "Hello!", $red);
		$bg->copy($graph, $x-50,$y-50,0,0, 100, 100);
	}


	# Save to a new file
	open my $out, '>:raw', "results/hsa${pImg}.png" or die "Cannot write output.png: $!";
	print $out $bg->png;
	close $out;
	<STDIN>;
}

sub load_metadata {
	my %args = (
		@_
	);
	my $file = $args{file};
	my $lipids = {};
	open(IN,$file);
	while(<IN>) {
		next if ($_ =~ /^Input/);
		my @l = split("\t",$_);
		next if (!$l[2]);
		my $lipid = {
			input_id => $l[0],
			kegg_id => $l[2],
			fa_type => $l[3],
			category => $l[4],
			main_class => $l[5],
			sub_class => $l[6]
		};
		$lipids->{data}->{$l[0]} = $lipid;
		$lipids->{metadata}-> {kegg2input} -> {$l[2]} -> {$l[0]} = {};
		$lipids->{metadata}-> {input2kegg} -> {$l[0]} -> {$l[2]} = {};
	}
	close(IN);
	return($lipids);
}

sub load_de_data {
	my %args = (
		@_
	);
	my $file = $args{file};
	my $metadata = $args{metadata};
	open(IN,$file);
	while(<IN>) {
		next if ($_ =~ /^Lipid/);
		my @l = split("\t",$_);
		next if (!$metadata->{data}->{$l[0]});
		$metadata->{data}->{$l[0]}->{dev} = $l[1];
		$metadata->{data}->{$l[0]}->{pval} = $l[5];
	}
	close(IN);
	return($metadata);
}

sub load_gmt {
	my %args = (
		@_
	);
	my $file = $args{file};
	my $metadata = $args{metadata};

	my $gmt;
	open(IN,$file);
	while(<IN>) {
		chomp;
		my ($pID,$pName,@metas) = split("\t",$_);
		$gmt->{metadata}->{id2name}->{$pID} = $pName;
		foreach my $meta (@metas) {
			if ($metadata->{metadata}->{kegg2input}->{$meta}) {
				$gmt->{pathway2meta}->{$pID}->{$meta} = {};
				$gmt->{meta2pathway}->{$meta}->{$pID} = {};
			}
		}
	}
	close(IN);
	return($gmt);
}

sub load_nodes {
	my %args = (
		@_
	);
	my $file = $args{file};
	my $data = $args{data};
	my @nodes;
	open(IN,$file) or die "Cannot open: $file\n";
	while(<IN>) {
		chomp;
		my @l = split("\t",$_);
		next if (!$data->{metadata}->{kegg2input}->{$l[3]});
		
		my ($x,$y) = split(",",$l[4]);
		my $node = {
			x => $x,
			y => $y,
			kegg_id => $l[3]
		};
		push(@nodes,$node);
	}
	close(IN);
	return(@nodes)
}

sub load_de_on_node {
	my %args = (
		@_
	);
	my $node = $args{node};
	my $data = $args{data};

	my $nodeID = $node->{kegg_id};
	my @nodeInputs = keys(%{$data->{metadata}->{kegg2input}->{$node->{kegg_id}}});

	foreach my $nodeInput (@nodeInputs) {
		if ($data->{data}->{$nodeInput}) {
			$node->{lipids}->{$nodeInput} = dclone $data->{data}->{$nodeInput}
		}
	}
}

sub make_complexes {
	my %args = (
		@_
	);
	my @nodes = @{$args{nodes}};
	my $index;
	my $n = 0;
	foreach my $node (@nodes) {
		my $coords = "$node->{x},$node->{y}";
		$index -> {$coords} -> {$n} = $node;
	}
	return($index);
}

sub sliced_crown_arc {
	my %args = (
		theta_start => 0,
		theta_end => 360,
		@_
	);
	my $graph = $args{graph};
	my $cx = $args{centerX};
	my $cy = $args{centerY};
	my $d1 = $args{diameter1};
	my $d2 = $args{diameter2};
	my @counts = @{$args{counts}};
	my @colors = @{$args{colors}};
	my $theta_start = $args{theta_start};
	my $theta_end = $args{theta_end};


	
	my $black = $graph->colorAllocate(0,0,0);

	my $total_counts = sum(@counts);
	my @xcents;
	foreach my $count (@counts) {
		push(@xcents,($count/$total_counts)*100);
	}
	if ($xcents[0] > $xcents[1]) {
		$xcents[0] = ceil($xcents[0]);
		$xcents[1] = floor($xcents[1]);
	} else {
		$xcents[0] = floor($xcents[0]);
		$xcents[1] = ceil($xcents[1]);
	}

	if ($xcents[0] < 20) {
		$xcents[0] = 20;
		$xcents[1] = 80;
	}
	if ($xcents[1] < 20) {
		$xcents[1] = 20;
		$xcents[0] = 80;
	}

	say("  up: ${xcents[0]}% - dn: ${xcents[1]}%")  if $debug;
	my $theta0 = min($theta_start,$theta_end);
	my $previous_theta = $theta0;
	my @thetas = ($theta0);
	my $invert;

	my $theta_slice = abs($theta_end - $theta_start);
	foreach my $xcent (@xcents) {
		my $theta = (($theta_slice * $xcent) / 100) + $previous_theta;
		$previous_theta = $theta;
		my $theta_r = ($theta/180)*3.14;
		push(@thetas,$theta_r)
	}

	$previous_theta = $theta0;
	my $n = 0;
	$graph->arc($cx,$cy,$d1,$d1,$theta_start,$theta_end,$black);
	$graph->arc($cx,$cy,$d2,$d2,$theta_start,$theta_end,$black);
	say Dumper \@thetas;
	foreach my $theta (@thetas) {
		my $p1y = $cy-($d1 / 2)* sin($theta);
		my $p1x = $cx-($d1 / 2)* cos($theta);
		my $p2x = $cx-($d2 / 2)* cos($theta);
		my $p2y = $cy-($d2 / 2)* sin($theta);
		say("  radiant: $theta") if $debug;
		say("  p1x: ${p1x} | p1y: ${p1y} | p2y: ${p2y}") if $debug;

		$graph->line($p1x, $p1y, $p2x, $p2y, $black);

		if ($theta > $theta0) {
			my $dm = $d1+(($d2-$d1)/2);
			my $theta_m = ($theta - $previous_theta)/2;
			my $pmx = $cx-($dm / 2)* cos($theta_m);
			my $pmy = $cy-($dm / 2)* sin($theta_m);
			if ($theta > 3.14) {
				$pmy = $cy+($cy-$pmy);
			}
			$graph->fill($pmx,$pmy,$colors[$n]);
			$n++;
		}
		$previous_theta = $theta;
	}
	
	return($graph);
}