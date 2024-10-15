use leptos::*;
use leptos_meta::*;
use leptos_router::*;

/// The main App component.
#[component]
pub fn App() -> impl IntoView {
    // Provides context for meta tags, stylesheets, etc.
    provide_meta_context();

    view! {
        // Injects a stylesheet into the document <head>
        // id=leptos means cargo-leptos will hot-reload this stylesheet
        <Stylesheet id="leptos" href="/pkg/server.css"/>
        // Sets the document title
        <Title text="DNCS"/>
        // Main content of the app
        <Router>
            <main>
                <Routes>
                    <Route path="" view=HomePage/>
                    <Route path="/*any" view=NotFound/>
                </Routes>
            </main>
        </Router>
    }
}

fn handle_form_submission(
    molecule: String,
    sequence: String,
    samples: i32,
    include_sidechain: bool,
) {
    // Handle form data
    println!("Molecule: {}", molecule);
    println!("Sequence: {}", sequence);
    println!("Samples: {}", samples);
    println!("Include Sidechain: {}", include_sidechain);
}

/// Renders the home page of your application.
#[component]
fn HomePage() -> impl IntoView {
    let (molecule, set_molecule) = create_signal(String::new());
    let (sequence, set_sequence) = create_signal(String::new());
    let (samples, set_samples) = create_signal(0);
    let (include_sidechain, set_include_sidechain) = create_signal(false);

    view! {
        <h1>"DNCS"</h1>
        <form class="config" on:submit=move |ev| {
            ev.prevent_default();
            let molecule_value = molecule.get();
            let sequence_value = sequence.get();
            let samples_value = samples.get();
            let include_sidechain_value = include_sidechain.get();

            handle_form_submission(
                molecule_value,
                sequence_value,
                samples_value,
                include_sidechain_value
            );
        }>
            <div class="parameter">
                <label for="molecule">
                    <b>"Molecule: "</b>
                </label>
                <input
                    type="text"
                    id="molecule"
                    name="molecule"
                    placeholder="Just Enter Some Name"
                    on:input=move |ev| {
                        set_molecule.set(event_target_value(&ev));
                    }
                />
            </div>
            <div class="parameter">
                <label for="sequence">
                    <b>"Sequence: "</b>
                </label>
                <input
                    type="text"
                    id="sequence"
                    name="sequence"
                    placeholder="Enter Amino Acid Sequence"
                    on:input=move |ev| {
                        set_sequence.set(event_target_value(&ev));
                    }
                />
            </div>
            <div class="parameter">
                <label for="samples">
                    <b>"Samples: "</b>
                </label>
                <input
                    type="number"
                    id="samples"
                    name="samples"
                    placeholder="No of Samples"
                    on:input=move |ev| {
                        if let Ok(value) = event_target_value(&ev).parse::<i32>() {
                            set_samples.set(value);
                        }
                    }
                />
            </div>
            <div class="parameter">
                <label for="include-sidechain">
                    <b>"Include Sidechain: "</b>
                </label>
                <input
                    type="checkbox"
                    id="include-sidechain"
                    name="include-sidechain"
                    on:change=move |ev| {
                        set_include_sidechain.set(event_target_checked(&ev));
                    }
                />
            </div>
            <div style="margin-top: 20px; text-align: center;">
                <button type="submit">"Generate"</button>
            </div>
        </form>
    }
}

/// 404 - Not Found page.
#[component]
fn NotFound() -> impl IntoView {
    #[cfg(feature = "ssr")]
    {
        // Set HTTP status code 404
        let resp = expect_context::<leptos_actix::ResponseOptions>();
        resp.set_status(actix_web::http::StatusCode::NOT_FOUND);
    }

    view! {
        <h1>"Not Found"</h1>
    }
}
