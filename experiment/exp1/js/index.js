///Don't touch lines 2-43.
var order = 1;
/// Helper function that shuffles an array. Don't touch.
var shuffle = function (array) {

	var currentIndex = array.length;
	var temporaryValue, randomIndex;
	while (0 !== currentIndex) {
		randomIndex = Math.floor(Math.random() * currentIndex);
		currentIndex -= 1;
		temporaryValue = array[currentIndex];
		array[currentIndex] = array[randomIndex];
		array[randomIndex] = temporaryValue;
	}

	return array;

};
///Helper functions to get random selection from array, and to remove elements from array. Don't touch.
function getRandomNonDuplicateSelection(arr, count, exclusionArr) {
  var selected = [];
  var available = arr.filter(item => !exclusionArr.includes(item));

  if (available.length < count) {
    throw new Error('Insufficient unique elements available for selection.');
  }

  for (var i = 0; i < count; i++) {
    var randomIndex = Math.floor(Math.random() * available.length);
    selected.push(available[randomIndex]);
    available.splice(randomIndex, 1);
  }
  return selected;
}

function removeFromArray(arr, elements) {
  for (var i = 0; i < elements.length; i++) {
    var index = arr.indexOf(elements[i]);
    if (index > -1) {
      arr.splice(index, 1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function generateStimulusList() {
  var numBlocks = 6;

  // Items 1..12
  var allItems = [];
  for (var i = 1; i <= 12; i++) allItems.push(i);

  // Fillers
  var grammaticalFillers = [];
  var ungrammaticalFillers = [];
  for (var j = 1; j <= 12; j++) {
    grammaticalFillers.push(8000 + j);    // 8001–8012
    ungrammaticalFillers.push(9000 + j);  // 9001–9012
  }

  // Randomly split items into 6 for cond1 and 6 for cond2
  shuffle(allItems);
  var cond1Items = allItems.slice(0, numBlocks);        // 6 items
  var cond2Items = allItems.slice(numBlocks, 2*numBlocks); // other 6 items

  // (Optional) reshuffle within each condition group so pairing is random
  shuffle(cond1Items);
  shuffle(cond2Items);

  var blocks = [];

  for (var b = 0; b < numBlocks; b++) {
    var itemCond1 = cond1Items[b];
    var itemCond2 = cond2Items[b];

    // Build critical IDs numerically: 1AB1 or 1AB2
    var crit1 = 1000 + (itemCond1 * 10) + 1; // condition 1
    var crit2 = 1000 + (itemCond2 * 10) + 2; // condition 2 (different item)

    // Pick one grammatical filler (no repeats)
    var gram = getRandomNonDuplicateSelection(grammaticalFillers, 1, [])[0];
    removeFromArray(grammaticalFillers, [gram]);

    // Pick one ungrammatical filler (no repeats)
    var ungram = getRandomNonDuplicateSelection(ungrammaticalFillers, 1, [])[0];
    removeFromArray(ungrammaticalFillers, [ungram]);

    var block = [crit1, crit2, gram, ungram];

    // Shuffle within block
    shuffle(block);

    blocks.push(block);
  }

  // Shuffle block order but keep blocks intact
  shuffle(blocks);
console.log(blocks);
  // Flatten
  var presentation_list = [];
  for (var k = 0; k < blocks.length; k++) {
    presentation_list = presentation_list.concat(blocks[k]);
  }

  return presentation_list;
}

// Store full list here (integers)
var presentation_list = generateStimulusList();
console.log(presentation_list);

////Now, we have a list of sentence IDs. We need to look up the actual sentences corresponding to the IDs in the stimuli.js file. 
const populated_list = presentation_list.map((integer) => all_stimuli.find((item) => item.unique_id === integer));
console.log(populated_list);




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////This section contain all the types of slides that will be used by the html file (e.g. practice slides, feedback slides, critial trial slides, demographic info slides, thank-you slide, etc.).////////
////////You can just modify the templates based on the need of your experiment. The most important thing is to update the "log_responses" functions, so that you log all the information you need.    ////////  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function make_slides(f) {
  var slides = {};  
  slides.i0 = slide({
     name : "i0",
     start: function() {
      exp.startT = Date.now();
     }
  });

  slides.instructions = slide({
    name : "instructions",
    button : function() {
      exp.go(); //use exp.go() if and only if there is no "present" data.
    }
  });

  slides.single_trial = slide({
    name: "single_trial",
    start: function() {
      $(".err").hide();
      $(".display_condition").html("You are in " + exp.condition + ".");
    },
    button : function() {
      response = $("#text_response").val();
      if (response.length == 0) {
        $(".err").show();
      } else {
        exp.data_trials.push({
          "trial_type" : "single_trial",
          "response" : response
        });
        exp.go(); //make sure this is at the *end*, after you log your data
      }
    },
  });

  slides.practice_slider = slide({
    name : "practice_slider",

    /* trial information for this block
     (the variable 'stim' will change between each of these values,
      and for each of these, present_handle will be run.) */
    present : [{"a": 1}],
    //this gets run only at the beginning of the block
    present_handle : function(stim) {
      $(".err").hide();
      $(".errgood").hide();
      this.stim = stim;
      $(".prompt").html("Who did you send the letter to?");
      this.init_sliders();
      exp.sliderPost = null; //erase current slider value
      exp.first_response_wrong = 0;
      exp.first_response_value = null;
      exp.attempts = 0;
    },
    button : function() {
      if (exp.sliderPost == null) {
        $(".err").show();
      } 
      else if (exp.sliderPost < 0.5) {
        exp.first_response_wrong = 1;
        exp.first_response_value =exp.sliderPost;
        exp.attempts = exp.attempts + 1;
        $(".errgood").show();
      }
      else {
        this.log_responses();
        /* use _stream.apply(this); if and only if there is
        "present" data. (and only *after* responses are logged) */
        _stream.apply(this);
      }
    },
    init_sliders : function() {
      utils.make_slider("#practice_slider_1", function(event, ui) {
        exp.sliderPost = ui.value;
      });
    },
    log_responses : function() {
      exp.data_trials.push({
        "response" : exp.sliderPost,
        "first_response_value": exp.first_response_value,
        "wrong_attempts": exp.attempts,
        "item_type" : "practice_good",
        "item_number": "practice_good",
        "trial_sequence_total": 0
      });

    }
  });


  slides.post_practice_1 = slide({
    name : "post_practice_1",
    button : function() {
      exp.go(); //use exp.go() if and only if there is no "present" data.
    }
  });


 

  slides.practice_slider_bad = slide({
    name : "practice_slider_bad",

    /* trial information for this block
     (the variable 'stim' will change between each of these values,
      and for each of these, present_handle will be run.) */
    present : [1],

  
    //this gets run only at the beginning of the block
    present_handle : function(stim) {
      $(".err").hide();
      $(".errbad").hide();
      $(".prompt").html("Who email did send he to?");
      this.init_sliders();
      exp.sliderPost = null; //erase current slider value
      exp.first_response_wrong = 0;
      exp.first_response_value = null;
      exp.attempts = 0;
    },
    button : function() {
      if (exp.sliderPost == null) {
        $(".err").show();
      } 
      else if (exp.sliderPost > 0.5) {
        exp.first_response_wrong = 1;
        exp.first_response_value = exp.sliderPost;
        exp.attempts = exp.attempts + 1;
        $(".errbad").show();
      }
      else {
        this.log_responses();
        /* use _stream.apply(this); if and only if there is
        "present" data. (and only *after* responses are logged) */
        _stream.apply(this);
      }
    },
    init_sliders : function() {
      utils.make_slider("#practice_slider_2", function(event, ui) {
        exp.sliderPost = ui.value;
        
      });
    },
    log_responses : function() {
      exp.data_trials.push({
        "response" : exp.sliderPost,
        "first_response_value": exp.first_response_value,
        "wrong_attempts": exp.attempts,
        "item_type" : "practice_bad",
        "item_number": "practice_bad",
        "trial_sequence_total": 0
      });

    }
  });

  slides.post_practice_2 = slide({
    name : "post_practice_2",
    button : function() {
      exp.go(); //use exp.go() if and only if there is no "present" data.
    }
  });


  slides.last_reminder = slide({
    name : "last_reminder",
    button : function() {
      exp.go(); //use exp.go() if and only if there is no "present" data.
    }
    
  });

 

  slides.one_slider = slide({
    name : "one_slider",

    /* trial information for this block
     (the variable 'stim' will change between each of these values,
      and for each of these, present_handle will be run.) */
    present : populated_list,
    
    //this gets run only at the beginning of the block
    present_handle : function(stim) {
      $(".err").hide();
      this.stim = stim; //I like to store this information in the slide so I can record it later.
      $(".target").html(stim.sentence);
      this.init_sliders()
      exp.sliderPost = null; //erase current slider value
    },

    button : function() {
      if (exp.sliderPost == null) {
        $(".err").show();
      } else {
        this.log_responses();

        /* use _stream.apply(this); if and only if there is
        "present" data. (and only *after* responses are logged) */
        _stream.apply(this);
      }
    },

    init_sliders : function() {
      utils.make_slider("#single_slider", function(event, ui) {
        exp.sliderPost = ui.value;
      });
    },

    log_responses : function() {
      exp.data_trials.push({
        // item-specific fields
        "response" : exp.sliderPost,
        "item": this.stim.item,
        "sentence": this.stim.sentence,
        "condition": this.stim.condition,
        "trial_sequence_total": order,
        "sentence_id": this.stim.unique_id
      });
      order = order + 1;
    }
  });

  slides.subj_info =  slide({
    name : "subj_info",
    submit : function(e){
      //if (e.preventDefault) e.preventDefault(); // I don't know what this means.
      exp.subj_data = {
        asses : $('input[name="assess"]:checked').val(),
        age : $("#age").val(),
        gender : $("#gender").val(),
        education : $("#education").val(),
        comments : $("#comments").val(),
        fairprice: $("#fairprice").val()
      };
      exp.go(); //use exp.go() if and only if there is no "present" data.
    }
  });

  slides.thanks = slide({
    name : "thanks",
    start : function() {
      exp.data= {
          "trials" : exp.data_trials,
          "catch_trials" : exp.catch_trials,
          "system" : exp.system,
          "subject_information" : exp.subj_data,
          "time_in_minutes" : (Date.now() - exp.startT)/60000
      };
      proliferate.submit(exp.data);
    }
  });

  return slides;
}

/// init ///
function init() {
  exp.trials = [];
  exp.catch_trials = [];
  //exp.condition = _.sample(["condition 1", "condition 2"]); //can randomize between subject conditions here
  exp.system = {
      Browser : BrowserDetect.browser,
      OS : BrowserDetect.OS,
      screenH: screen.height,
      screenUH: exp.height,
      screenW: screen.width,
      screenUW: exp.width
    };
  //blocks of the experiment:
  exp.structure=["i0", "instructions", "practice_slider", "post_practice_1", "practice_slider_bad", "post_practice_2", "last_reminder", 'one_slider', 'subj_info', 'thanks'];

  exp.data_trials = [];
  //make corresponding slides:
  exp.slides = make_slides(exp);

  exp.nQs = utils.get_exp_length(); //this does not work if there are stacks of stims (but does work for an experiment with this structure)
                    //relies on structure and slides being defined

  $('.slide').hide(); //hide everything

  //make sure turkers have accepted HIT (or you're not in mturk)
  $("#start_button").click(function() {
    exp.go();
  });

  exp.go(); //show first slide
}
