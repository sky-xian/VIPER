/**--------------------------------------
@author: Mahesh Vangala
@email: vangalamaheshh@gmail.com
@data: June, 30, 2016 
----------------------------------------*/

// vim: syntax=javascript tabstop=4 expandtab

var data = [{
		name: 'alignment',
		url: './partials/alignment.html',
		children: [{
			name: 'genome',
            url: './partials/genome_align.html'
		},{
			name: 'rRNA',
            url: './partials/rRNA_align.html'
		}]
	},{
		name: 'readQC',
        url: './partials/readQC.html',
        children: [{
            name: 'read_distribution',
            url: './partials/read_distrib.html'
        }]
	},{
		name: 'DE',
        children: [{
            name: 'volcano_plots'
        }]
	}];

$(document).ready(function(){
	$(function(){
        var $tree = $('#m-left');
		$tree.tree({ 
			data: data,
			autoOpen: true,
			dragAndDrop: false
		}).bind('tree.click', function(event){
			var node = event.node;
			var theURL = node.url;
			if( theURL ){ $('#m-right').load(theURL); }
		});
        var node = $tree.tree('getNodeByName', 'alignment');
        $tree.tree('selectNode',node); 
        $('#m-right').load(node.url)
	});
});
